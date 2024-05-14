clear
dim = 3; % State dimension
numElement = 2*dim;
numMea = 200;   % Number of measurements
numVair = 2; 
monTimes = 2e3;
numPoc = 1;

numGau = 3;
varPhiInit = zeros(numGau,1);

numSampTol = [20 50 100 200 500 800 1000];

lengthSamp = length(numSampTol);

devSampCKF = zeros(2,lengthSamp);
devSampSCKF = zeros(2,lengthSamp);
devSampPF = zeros(2,lengthSamp);
crbSamp = zeros(2,lengthSamp);
crbSampSCKF = zeros(2,lengthSamp);

timDelta = 0.05;    % time interval

delPosAP = 100; %Sensor interval
squSize = 2e3; %C_V2X Network Size

% variance of observations
varDis = 1e-2; % Variance of distance measurements
varSpe = 1; % Variance of radial velocity measurements




devLocUpdateFinal = zeros(numVair,numMea);
devLocUpdateMag_Final = zeros(numVair,numMea);
devLocKal_Final = zeros(numVair,numMea);
devLocPF_Final = zeros(numVair,numMea);
devLocSCKF_Final = zeros(numVair,numMea);
devLocUKF_Final = zeros(numVair,numMea);

% state transition matrix
statTransSca = [1 0;timDelta 1];
statTrans = kron(statTransSca,eye(dim));
statInSca = [timDelta 1/2*timDelta^2]';
statIn = kron(statInSca,eye(dim));

% Acceleration settings
varPhiInit(1) = 1e-1;
varPhiInit(2) = 1e-1;
varPhiInit(3) = 0;

phiMeanInit1 = [2 2 0]';
phiMeanInit2 = [-2 -2 0]';
phiMeanInit3 = zeros(dim,1);

Prob1 = 0.2;
Prob2 = 0.2;
Prob3 = 1 - Prob1 - Prob2;
Prob = [Prob1 Prob2 Prob3];

acc = zeros(dim,numMea+1);
accMean = zeros(3*dim,numMea+1);
varAcc1 = zeros(dim,numMea+1);
varAcc2 = zeros(dim,numMea+1);
varAcc3 = zeros(dim,numMea+1);
varAcc1(:,1) = varPhiInit(1)*ones(dim,1);
varAcc2(:,1) = varPhiInit(2)*ones(dim,1);
varAcc3(:,1) = varPhiInit(3)*ones(dim,1);

phiPdf1 = phiMeanInit1 + sqrt(varPhiInit(1))*randn(dim,1);
phiPdf2 = phiMeanInit2 + sqrt(varPhiInit(2))*randn(dim,1);
phiPdf3 = phiMeanInit3 + sqrt(varPhiInit(3))*randn(dim,1);

chi = randsrc(dim,1,[1,2,3;Prob1,Prob2,Prob3]);
phiPdf = zeros(dim,numMea+1);

phiPdf(:,1) = phiPdf1.*(chi==1)+phiPdf2.*(chi==2)+phiPdf3.*(chi==3);
acc(:,1) = phiPdf(:,1);
accMean(:,1) = [phiMeanInit1' phiMeanInit2' phiMeanInit3']';

varStatFir = varPhiInit(1);
nodeLocUpdate = zeros(numElement,numMea+1);
nodeLocUpdate(:,1) = zeros(numElement,1);


crbFinal = zeros(2,numMea);
crbIniIdx = [0.04 0.01]';
crbIni = kron(diag(crbIniIdx),eye(dim)); 

crbFinal(:,1) = crbIniIdx;


crbFinalSCKF = zeros(2,numMea);
crbFinalSCKF(:,1) = crbIniIdx;

varEKFFinal = zeros(2,numMea);
varEKFFinal(:,1) = crbIniIdx;

nodeCKFEst(:,1) = nodeLocUpdate(:,1) + sqrtm(crbIni)*randn(numElement,1);

nodeLocEst = zeros(dim,numMea+1);
nodeLocEst(:,1) = nodeCKFEst(dim+1:end,1);

nodeLocEst1 = zeros(dim,numMea+1);
nodeLocEst1(:,1) = nodeLocUpdate(dim+1:2*dim,1);

devLocUpdateFinal(:,1) = crbIniIdx;
devLocUpdateMag_Final(:,1) = crbIniIdx;
devLocKal_Final(:,1) = crbIniIdx;
devLocSCKF_Final(:,1) = crbIniIdx;


for sampIdx = 1:lengthSamp
    for monIdx = 1:1:monTimes
        sampIdx
        monIdx
        % genie path
        for meaIdx = 2:1:numMea+1
            chi = randsrc(dim,1,[1,2,3;Prob1,Prob2,Prob3]);
            phiPdf1 = phiMeanInit1 + sqrt(varPhiInit(1))*randn(dim,1);
            phiPdf2 = phiMeanInit2 + sqrt(varPhiInit(2))*randn(dim,1);
            phiPdf3 = phiMeanInit3 + sqrt(varPhiInit(3))*randn(dim,1);
            accMean(:,meaIdx) = accMean(:,1);
            acc(:,meaIdx) = phiPdf1.*(chi==1)+phiPdf2.*(chi==2)+phiPdf3.*(chi==3);

            varAcc1(:,meaIdx) = varAcc1(:,1);
            varAcc2(:,meaIdx) = varAcc2(:,1);
            varAcc3(:,meaIdx) = varAcc3(:,1);

            nodeLocUpdate(:,meaIdx) = statTrans*nodeLocUpdate(:,meaIdx-1) + statIn*acc(:,meaIdx-1);
        end

        % APs' position
        [posAP,sign,numAP] = CV2X_Network(nodeLocUpdate,dim,squSize,delPosAP);

        monTimesReal = monTimes;

        %% Initialization
        flag = 0;
        tic
        % CRB-KF - initialize the variables
        crbTol = zeros(numElement,numMea);
        crbTol1 = zeros(numElement,numMea);

        devLocUpdate = zeros(2,numMea);
        nodeLoc_last = nodeCKFEst(:,1);
        R_tol = zeros(2*dim,2*dim,numMea+1);
        M_tol = zeros(2*dim,2*dim,numMea+1);
        S_tol = zeros(2*dim,2*dim,numMea+1);
        R_tol(:,:,1) = crbIni;
        M_tol(:,:,1) = crbIni;
        S_tol(:,:,1) = crbIni;
        R_tol(:,:,2) = statTrans*crbIni*statTrans' + statIn*diag(varAcc1(:,1))*statIn';
        M_tol(:,:,2) = statTrans*crbIni*statTrans' + statIn*diag(varAcc2(:,1))*statIn';
        S_tol(:,:,2) = statTrans*crbIni*statTrans' + statIn*diag(varAcc3(:,1))*statIn';


        %   SCKF - initialize the variables
        statForeSCKF = zeros(2*dim,numMea);
        statForeSCKF(:,1) = randn(2*dim,1);
        statEstSCKF = zeros(2*dim,numMea+1);
        statEstSCKF(:,1) = nodeCKFEst(:,1);
        devLocSCKF = zeros(2,numMea);
        R_tol_SCKF = zeros(2*dim,2*dim,numMea+1);
        M_tol_SCKF = zeros(2*dim,2*dim,numMea+1);
        S_tol_SCKF = zeros(2*dim,2*dim,numMea+1);
        R_tol_SCKF(:,:,1) = crbIni;
        M_tol_SCKF(:,:,1) = crbIni;
        S_tol_SCKF(:,:,1) = crbIni;

        R_tol_SCKF(:,:,2) = statTrans*crbIni*statTrans' + statIn*diag(varAcc1(:,1))*statIn';
        M_tol_SCKF(:,:,2) = statTrans*crbIni*statTrans' + statIn*diag(varAcc2(:,1))*statIn';
        S_tol_SCKF(:,:,2) = statTrans*crbIni*statTrans' + statIn*diag(varAcc3(:,1))*statIn';

        P_SCKF = Prob1*R_tol_SCKF(:,:,1) + Prob2*M_tol_SCKF(:,:,1) + Prob3*S_tol_SCKF(:,:,1);

        %   PF - initialize the variables
        numSamplePF = numSampTol(sampIdx);
        nodeLocPF_Sample = zeros(numElement,numSamplePF);
        for PFstep = 1:1:numSamplePF
            nodeLocPF_Sample(:,PFstep) = nodeLocUpdate(:,1) + sqrt(varStatFir)*randn(numElement,1);
        end
        meaUpdatePF = zeros(2,numSamplePF);
        preWeightPF = zeros(1,numSamplePF);
        nodeEstPF = zeros(numElement,numMea);
        devLocPF = zeros(2,numMea);


        relDisTol = zeros(numMea,2);
        posAncTol = zeros(dim,numMea);

        for meaIdx = 2:1:numMea

            % Acceleration Initialization
            R_acc = diag(varAcc1(:,meaIdx));
            M_acc = diag(varAcc2(:,meaIdx));
            S_acc = diag(varAcc3(:,meaIdx));

            delAccMean = Prob1*accMean(1:dim,meaIdx) + Prob2*accMean(dim+1:2*dim,meaIdx) + Prob3*accMean(2*dim+1:3*dim,meaIdx);
            delAccVar = Prob1*(accMean(1:dim,meaIdx)*accMean(1:dim,meaIdx)'+ R_acc) +  ...
                Prob2*(accMean(dim+1:2*dim,meaIdx)*accMean(dim+1:2*dim,meaIdx)' + M_acc) + ...
                Prob3*(accMean(2*dim+1:3*dim,meaIdx)*accMean(2*dim+1:3*dim,meaIdx)' +S_acc) - ...
                delAccMean*delAccMean';

            % Randomly find one of the k closest APs
            relDisIdx = zeros(numAP,1);
            for idx = 1:numAP
                relDisIdx(idx) = norm(nodeLocUpdate(dim+1:end,meaIdx)-posAP(:,idx));
            end

            [~,sign] = mink(relDisIdx,10);
            sign = randsrc(1,1,sign');

            %% Measurements
            posAnc = posAP(:,sign); % AP position
            posAncTol(:,meaIdx) = posAnc;
            varMeas = [varDis;varSpe];
            noiMeas = sqrt(varMeas).*randn(2,1);
            measTol = mea_trans(nodeLocUpdate(:,meaIdx),posAnc,dim); % Distance and radial speed (measurement)
            relDisTol(meaIdx,1) = measTol(1);
            relDisTol(meaIdx,2) = sign;
            measTol = measTol + noiMeas;

            %% Estimation

            %PF
            accMean1 = accMean(1:dim,meaIdx);
            accMean2 = accMean(dim+1:2*dim,meaIdx);
            accMean3 = accMean(2*dim+1:3*dim,meaIdx);
            for PFstep = 1:1:numSamplePF
                chi = randsrc(dim,1,[1,2,3;Prob1,Prob2,Prob3]);
                pdf1 = accMean1 + sqrt(varAcc1(:,meaIdx)).*randn(dim,1);
                pdf2 = accMean2 + sqrt(varAcc2(:,meaIdx)).*randn(dim,1);
                pdf3 = accMean3 + sqrt(varAcc3(:,meaIdx)).*randn(dim,1);
                accPF = pdf1.*(chi==1)+pdf2.*(chi==2)+pdf3.*(chi==3);
                nodeLocPF_Sample(:,PFstep) = statTrans*nodeLocPF_Sample(:,PFstep) + statIn*accPF;
                meaUpdatePF(:,PFstep) = mea_trans(nodeLocPF_Sample(:,PFstep),posAnc,dim);
                Q = varMeas.*eye(2);
                preWeightPF(PFstep) = 1e-9 + 1/sqrt((2*pi)^2*det(Q)) * exp(-1/2*(measTol - meaUpdatePF(:,PFstep))'*(Q\(measTol - meaUpdatePF(:,PFstep))));
            end
            preWeightPF = preWeightPF/sum(preWeightPF);
            for PFstep = 1:1:numSamplePF
                nodeLocPF_Sample(:,PFstep) = nodeLocPF_Sample(:,find(rand <= cumsum(preWeightPF,2),1));   % 粒子权重大的将多得到后代
            end

            nodeEstPF(:,meaIdx) = sum(nodeLocPF_Sample,2)/numSamplePF;
            devLocPF(1,meaIdx) = norm(nodeEstPF(1:dim,meaIdx) - nodeLocUpdate(1:dim,meaIdx))^2;
            devLocPF(2,meaIdx) = norm(nodeEstPF(dim+1:2*dim,meaIdx) - nodeLocUpdate(dim+1:2*dim,meaIdx))^2;


            %SCKF
            statForeSCKF(:,meaIdx) = statTrans*statEstSCKF(:,meaIdx-1) + statIn*delAccMean;
            P_fore_SCKF = statTrans*P_SCKF*statTrans'+ statIn*delAccVar*statIn';
            h = partial_h(statForeSCKF(:,meaIdx),posAnc,dim);
            K = P_fore_SCKF*h'/(h*P_fore_SCKF*h'+ kron(diag(varMeas),eye(numPoc)));
            statEstSCKF(:,meaIdx) = statForeSCKF(:,meaIdx) + K*(measTol - mea_trans(statForeSCKF(:,meaIdx),posAnc,dim));
            [P_SCKF,~] = CRBSamp(statEstSCKF(:,meaIdx),posAnc,varMeas,statEstSCKF(:,meaIdx-1),statTrans,statIn,accMean(:,meaIdx),R_tol_SCKF(:,:,meaIdx),M_tol_SCKF(:,:,meaIdx),S_tol_SCKF(:,:,meaIdx),Prob,dim,numSampTol(sampIdx));
            crbSCKFIdx = diag(P_SCKF);
            crbFinalSCKF(1,meaIdx) = crbFinalSCKF(1,meaIdx) + sum(crbSCKFIdx(1:dim));
            crbFinalSCKF(2,meaIdx) = crbFinalSCKF(2,meaIdx) + sum(crbSCKFIdx(dim+1:2*dim));

            devLocSCKF(1,meaIdx) = norm(statEstSCKF(1:dim,meaIdx) - nodeLocUpdate(1:dim,meaIdx))^2;
            devLocSCKF(2,meaIdx) = norm(statEstSCKF(dim+1:2*dim,meaIdx) - nodeLocUpdate(dim+1:2*dim,meaIdx))^2;
            R_tol_SCKF(:,:,meaIdx+1) = statTrans*P_SCKF*statTrans' + statIn*R_acc*statIn';
            M_tol_SCKF(:,:,meaIdx+1) = statTrans*P_SCKF*statTrans' + statIn*M_acc*statIn';
            S_tol_SCKF(:,:,meaIdx+1) = statTrans*P_SCKF*statTrans' + statIn*S_acc*statIn';



            %CRB_KF
            nodeLoc_pre = nodeLoc_last;
            [nodeLoc_pre,crb] = CRB_KF_Samp(measTol,nodeLoc_pre,posAnc,dim,varMeas,nodeLoc_last,R_tol,M_tol,S_tol,statTrans,statIn,accMean,Prob,meaIdx,numSampTol(sampIdx));
            locErr(1) = norm(nodeLoc_pre(1:dim) - nodeLocUpdate(1:dim,meaIdx))^2;
            locErr(2) = norm(nodeLoc_pre(dim+1:2*dim) - nodeLocUpdate(dim+1:2*dim,meaIdx))^2;
            devLocUpdate(:,meaIdx) = locErr;
            crbTol(:,meaIdx) = diag(crb);
            crbFinal(1,meaIdx) = crbFinal(1,meaIdx) + sum(crbTol(1:dim,meaIdx));
            crbFinal(2,meaIdx) = crbFinal(2,meaIdx) + sum(crbTol(dim+1:2*dim,meaIdx));
            R_tol(:,:,meaIdx+1) = statTrans*crb*statTrans' + statIn*R_acc*statIn';
            M_tol(:,:,meaIdx+1) = statTrans*crb*statTrans' + statIn*M_acc*statIn';
            S_tol(:,:,meaIdx+1) = statTrans*crb*statTrans' + statIn*S_acc*statIn';
            nodeLoc_last = nodeLoc_pre;
            nodeCKFEst(:,meaIdx) = nodeLoc_pre;
            nodeLocEst(:,meaIdx) = nodeLoc_pre(dim+1:2*dim);



        end
        devLocUpdateFinal = devLocUpdateFinal + devLocUpdate;
        devLocPF_Final = devLocPF_Final + devLocPF;
        devLocSCKF_Final = devLocSCKF_Final + devLocSCKF;
%         fprintf('the numer of Monte Carlo is%6.1f\n',monIdx)
        toc
    end

    devLocPF_Final = sqrt(devLocPF_Final/monTimesReal);
    devLocUpdateFinal = sqrt(devLocUpdateFinal/monTimesReal);
    devLocSCKF_Final = sqrt(devLocSCKF_Final/monTimesReal);
    crbFinal = sqrt(crbFinal/monTimesReal);
    crbFinalSCKF = sqrt(crbFinalSCKF/monTimesReal);

    devSampCKF(:,sampIdx) = devLocUpdateFinal(:,numMea);
    devSampSCKF(:,sampIdx) = devLocSCKF_Final(:,numMea);
    devSampPF(:,sampIdx) = devLocPF_Final(:,numMea);
    crbSamp(:,sampIdx) = crbFinal(:,numMea);
    crbSampSCKF(:,sampIdx) = crbFinalSCKF(:,numMea);
    

end



figure
set(gcf,'position',[796,130,786,741])
subplot('Position',[0.085,0.587594128333916,0.88,0.35])
plot(numSampTol,devSampCKF(1,:),'-ob','linewidth',2.3);
hold on
plot(numSampTol,devSampSCKF(1,:),'-or','linewidth',2.3);
hold on
plot(numSampTol,crbSamp(1,:),'r','linestyle','--','marker','o','linewidth',2.3);
hold on
plot(numSampTol,crbSampSCKF(1,:),'r','linestyle','--','marker','o','linewidth',2.3);
hold on
plot(numSampTol,devSampPF(1,:),'color',[0.93,0.69,0.13],'marker','h','linewidth',2.3);
hold on
ylabel('RMSE (m/s)')
xlabel('Sample size')
legend('CRB-KF','SCKF','CRB','CRB-SCKF','PF')
title('Velocity Estimation')
set(gca,'FontSize',13);
grid on
hold on

subplot('Position',[0.085,0.091644204851752,0.88,0.35])
plot(numSampTol,devSampCKF(2,:),'-ob','linewidth',2.3);
hold on
plot(numSampTol,devSampSCKF(2,:),'-or','linewidth',2.3);
hold on
plot(numSampTol,crbSamp(2,:),'r','linestyle','--','marker','o','linewidth',2.3);
hold on
plot(numSampTol,crbSampSCKF(2,:),'r','linestyle','--','marker','o','linewidth',2.3);
hold on
plot(numSampTol,devSampPF(2,:),'color',[0.93,0.69,0.13],'marker','h','linewidth',2.3);
hold on
ylabel('RMSE (m)')
xlabel('Sample size')
legend('CRB-KF','SCKF','CRB','CRB-SCKF','PF')
title('Position Estimation')
set(gca,'FontSize',13);
grid on
hold on
name = ['Delta_Samp;MonTimes=' num2str(monTimesReal) '.fig'];
saveas(gcf,name);