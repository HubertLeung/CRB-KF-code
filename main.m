clear

dim = 3; % State dimension
numElement = 2*dim;
numMea = 500;   % Number of measurements
numVair = 2;    %
monTimes = 1;5e3;
numPoc = 1;

numGau = 3;
varPhiInit = zeros(numGau,1);

timDelta = 0.05;    % time interval

delPosAP = 100; %Sensor interval
squSize = 2e3; %C_V2X Network Size

% Variance of observations
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

nodeLocEst_CubKF = zeros(dim,numMea+1);
nodeLocEst_CubKF(:,1) = nodeLocUpdate(dim+1:2*dim,1);

nodeLocEst_SR = zeros(dim,numMea+1);
nodeLocEst_SR(:,1) = nodeLocUpdate(dim+1:2*dim,1);

devLocUpdateFinal(:,1) = crbIniIdx;
devLocCubKF_Final(:,1) = crbIniIdx;
devLocKal_Final(:,1) = crbIniIdx;
devLocSCKF_Final(:,1) = crbIniIdx;
devLocSRCKF_Final(:,1) = crbIniIdx;

for monIdx = 1:1:monTimes
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

    % Cub_KF --- Initialize
    devLocCubKF = zeros(2,numMea);
    nodeLoc_last_CubKF = nodeCKFEst(:,1);
    P_CubKF = Prob1*R_tol(:,:,1) + Prob2*M_tol(:,:,1) + Prob3*S_tol(:,:,1);


    % SRCKF --- Initialize
    devLocSR = zeros(2,numMea);
    nodeLoc_last_SR = nodeCKFEst(:,1);
    PSR = Prob1*R_tol(:,:,1) + Prob2*M_tol(:,:,1) + Prob3*S_tol(:,:,1);
    root_PSR = chol(PSR,'lower');

    %   EKF - initialize the variables
    statForeEKF = zeros(2*dim,numMea);
    statForeEKF(:,1) = randn(2*dim,1);
    statEstEKF = zeros(2*dim,numMea+1);
    statEstEKF(:,1) = nodeCKFEst(:,1);
    devLocKal = zeros(2,numMea);
    R_tol_Kal = crbIni;
    M_tol_Kal = crbIni;
    S_tol_Kal = crbIni;
    P = Prob1*R_tol_Kal + Prob2*M_tol_Kal + Prob3*S_tol_Kal;

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

    %   UKF - initialize the variables
    statForeUKF = zeros(2*dim,numMea);
    statForeUKF(:,1) = randn(2*dim,1);
    statEstUKF = zeros(2*dim,numMea+1);
    statEstUKF(:,1) = nodeCKFEst(:,1);
    devLocUKF = zeros(2,numMea);
    R_tol_UKF = statTrans*crbIni*statTrans' + statIn*diag(varAcc1(:,1))*statIn';
    M_tol_UKF = statTrans*crbIni*statTrans' + statIn*diag(varAcc2(:,1))*statIn';
    S_tol_UKF = statTrans*crbIni*statTrans' + statIn*diag(varAcc3(:,1))*statIn';
    PUKF = Prob1*R_tol_UKF + Prob2*M_tol_UKF + Prob3*S_tol_UKF;

    relDisTol = zeros(numMea,2);
    posAncTol = zeros(dim,numMea);

    for meaIdx = 2:1:numMea

        % % Acceleration Initialization
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
        posAnc = posAP(:,sign); %% Position of the AP
        posAncTol(:,meaIdx) = posAnc;
        varMeas = [varDis;varSpe];
        noiMeas = sqrt(varMeas).*randn(2,1);
        measTol = mea_trans(nodeLocUpdate(:,meaIdx),posAnc,dim); % Distance and radial speed (measurement)
        relDisTol(meaIdx,1) = measTol(1);
        relDisTol(meaIdx,2) = sign;
        measTol = measTol + noiMeas;

        %EKF
        statForeEKF(:,meaIdx) = statTrans*statEstEKF(:,meaIdx-1) + statIn*delAccMean;
        P_fore = statTrans*P*statTrans'+ statIn*delAccVar*statIn';
        h = partial_h(statForeEKF(:,meaIdx),posAnc,dim);
        K = P_fore*h'/(h*P_fore*h'+ kron(diag(varMeas),eye(numPoc)));
        statEstEKF(:,meaIdx) = statForeEKF(:,meaIdx) + K*(measTol - mea_trans(statForeEKF(:,meaIdx),posAnc,dim));
        P = (eye(2*dim)-K*h)*P_fore;
        varEKFIdx = diag(P);
        varEKFFinal(1,meaIdx) = varEKFFinal(1,meaIdx) + sum(varEKFIdx(1:dim));
        varEKFFinal(2,meaIdx) = varEKFFinal(2,meaIdx) + sum(varEKFIdx(dim+1:2*dim));

        devLocKal(1,meaIdx) = norm(statEstEKF(1:dim,meaIdx) - nodeLocUpdate(1:dim,meaIdx))^2;
        devLocKal(2,meaIdx) = norm(statEstEKF(dim+1:2*dim,meaIdx) - nodeLocUpdate(dim+1:2*dim,meaIdx))^2;



        %SCKF
        statForeSCKF(:,meaIdx) = statTrans*statEstSCKF(:,meaIdx-1) + statIn*delAccMean;
        P_fore_SCKF = statTrans*P_SCKF*statTrans'+ statIn*delAccVar*statIn';
        h = partial_h(statForeSCKF(:,meaIdx),posAnc,dim);
        K = P_fore_SCKF*h'/(h*P_fore_SCKF*h'+ kron(diag(varMeas),eye(numPoc)));
        statEstSCKF(:,meaIdx) = statForeSCKF(:,meaIdx) + K*(measTol - mea_trans(statForeSCKF(:,meaIdx),posAnc,dim));
        [P_SCKF,~] = CRB(statEstSCKF(:,meaIdx),posAnc,varMeas,statEstSCKF(:,meaIdx-1),statTrans,statIn,accMean(:,meaIdx),R_tol_SCKF(:,:,meaIdx),M_tol_SCKF(:,:,meaIdx),S_tol_SCKF(:,:,meaIdx),Prob,dim);
        crbSCKFIdx = diag(P_SCKF);
        crbFinalSCKF(1,meaIdx) = crbFinalSCKF(1,meaIdx) + sum(crbSCKFIdx(1:dim));
        crbFinalSCKF(2,meaIdx) = crbFinalSCKF(2,meaIdx) + sum(crbSCKFIdx(dim+1:2*dim));

        devLocSCKF(1,meaIdx) = norm(statEstSCKF(1:dim,meaIdx) - nodeLocUpdate(1:dim,meaIdx))^2;
        devLocSCKF(2,meaIdx) = norm(statEstSCKF(dim+1:2*dim,meaIdx) - nodeLocUpdate(dim+1:2*dim,meaIdx))^2;
        R_tol_SCKF(:,:,meaIdx+1) = statTrans*P_SCKF*statTrans' + statIn*R_acc*statIn';
        M_tol_SCKF(:,:,meaIdx+1) = statTrans*P_SCKF*statTrans' + statIn*M_acc*statIn';
        S_tol_SCKF(:,:,meaIdx+1) = statTrans*P_SCKF*statTrans' + statIn*S_acc*statIn';


        % UKF
        %         % Parameter settings
        %         a = 1e-2; %(0,1]
        %         b = 2; % =2
        %         k = 1; % >=0
        %
        %         weigSumMean_1 = Prob1*accMean(1:dim,meaIdx) + Prob2*accMean(dim+1:2*dim,meaIdx) + Prob3*accMean(2*dim+1:3*dim,meaIdx);
        %         DeltaVar_1 = Prob1*(accMean(1:dim,meaIdx)*accMean(1:dim,meaIdx)'+ R_acc) +  ...
        %             Prob2*(accMean(dim+1:2*dim,meaIdx)*accMean(dim+1:2*dim,meaIdx)' + M_acc) + ...
        %             Prob3*(accMean(2*dim+1:3*dim,meaIdx)*accMean(2*dim+1:3*dim,meaIdx)' +S_acc) - ...
        %             weigSumMean_1*weigSumMean_1';
        %
        %         QUKF = statIn*DeltaVar_1*statIn';%Prob1*statIn*R_acc*statIn' + Prob2*statIn*M_acc*statIn' + Prob3*statIn*S_acc*statIn';
        %         RUKF = diag(varMeas);
        %         [statEstUKF(:,meaIdx), PUKF] = UKFStep(statEstUKF(:,meaIdx-1), measTol, PUKF, QUKF, RUKF, a, k, b, dim,Prob,accMean(:,meaIdx),statTrans,statIn,posAnc,numVair);
        % %         [statEstUKF(:,meaIdx), PUKF] = UKFStep1(statEstUKF(:,meaIdx-1), measTol, PUKF, QUKF, RUKF, a, k, b, dim,Prob,accMean(:,meaIdx),statTrans,statIn,posAnc,numVair,numPoc);

%         [statEstUKF(:,meaIdx), PUKF] = UKF_weight(varMeas,dim,posAnc,measTol,delAccMean,delAccVar,statTrans,statIn,statEstUKF(:,meaIdx-1),PUKF,numElement);
%         devLocUKF(1,meaIdx) = norm(statEstUKF(1:dim,meaIdx) - nodeLocUpdate(1:dim,meaIdx))^2;
%         devLocUKF(2,meaIdx) = norm(statEstUKF(dim+1:2*dim,meaIdx) - nodeLocUpdate(dim+1:2*dim,meaIdx))^2;


        % CRB_KF
        nodeLoc_pre = nodeLoc_last;
        [nodeLoc_pre,crb] = CRB_KF(measTol,nodeLoc_pre,posAnc,dim,varMeas,nodeLoc_last,R_tol,M_tol,S_tol,statTrans,statIn,accMean,Prob,meaIdx); 
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

        % Cub_KF
        [stateUpd, P_CubKF] = cub_KF(varMeas,dim,posAnc,measTol,delAccMean,delAccVar,statTrans,statIn,nodeLoc_last_CubKF,P_CubKF,numElement);
        nodeLoc_last_CubKF = stateUpd;
        nodeLocEst_CubKF(:,meaIdx) = stateUpd(dim+1:2*dim);

        devLocCubKF(1,meaIdx) = norm(stateUpd(1:dim) - nodeLocUpdate(1:dim,meaIdx))^2;
        devLocCubKF(2,meaIdx) = norm(stateUpd(dim+1:2*dim) - nodeLocUpdate(dim+1:2*dim,meaIdx))^2;

        % SRCKF
        [state_SRCKF, PSR, root_PSR] = SRCKF(varMeas,dim,posAnc,measTol,delAccMean,delAccVar,statTrans,statIn,nodeLoc_last_SR,root_PSR,numElement);
        nodeLoc_last_SR = state_SRCKF;
        nodeLocEst_SR(:,meaIdx) = state_SRCKF(dim+1:2*dim);

        devLocSR(1,meaIdx) = norm(state_SRCKF(1:dim) - nodeLocUpdate(1:dim,meaIdx))^2;
        devLocSR(2,meaIdx) = norm(state_SRCKF(dim+1:2*dim) - nodeLocUpdate(dim+1:2*dim,meaIdx))^2;


    end

    devLocUpdateFinal = devLocUpdateFinal + devLocUpdate;
    devLocCubKF_Final = devLocCubKF_Final + devLocCubKF;
    devLocKal_Final = devLocKal_Final + devLocKal;
    devLocUKF_Final = devLocUKF_Final + devLocUKF;
    devLocSCKF_Final = devLocSCKF_Final + devLocSCKF;
    devLocSRCKF_Final = devLocSRCKF_Final + devLocSR;
    fprintf('the numer of Monte Carlo is%6.1f\n',monIdx)

    %     end
    toc
end

devLocUpdateFinal = sqrt(devLocUpdateFinal/monTimesReal);
devLocSCKF_Final = sqrt(devLocSCKF_Final/monTimesReal);
devLocCubKF_Final = sqrt(devLocCubKF_Final/monTimesReal);
crbFinal = sqrt(crbFinal/monTimesReal);
crbFinalSCKF = sqrt(crbFinalSCKF/monTimesReal);
varEKFFinal = sqrt(varEKFFinal/monTimesReal);
devLocKal_Final = sqrt(devLocKal_Final/monTimesReal);
devLocUKF_Final = sqrt(devLocUKF_Final/monTimesReal);
devLocSRCKF_Final = sqrt(devLocSRCKF_Final/monTimesReal);

% LL = length(posAncTol(1,:));
% figure
% for idx = 1:LL
%     plot(posAncTol(1,idx),posAncTol(2,idx),'o');
%     hold on
% end
%
% figure
% for idx = 1:numAP
%     plot(posAP(1,idx),posAncTol(2,idx),'o');
%     hold on
% end

Time = (2*timDelta:timDelta:numMea*timDelta);

figure
set(gcf,'position',[796,130,786,741])
subplot('Position',[0.085,0.587594128333916,0.88,0.35])
plot(Time,devLocUpdateFinal(1,2:end),'b','linewidth',2.3);
hold on
plot(Time,devLocSCKF_Final(1,2:end),'color',[0.69,0.19,0.38],'linewidth',2.3);
hold on
plot(Time,devLocSRCKF_Final(1,2:end),'m','linewidth',4);
hold on
plot(Time,devLocCubKF_Final(1,2:end),'g','linewidth',3);
hold on
plot(Time,devLocKal_Final(1,2:end),'k','linestyle','-.','linewidth',2.3);
hold on
plot(Time,real(crbFinal(1,2:end)),'r','linestyle','--','linewidth',2.3);
hold on
plot(Time,real(crbFinalSCKF(1,2:end)),'color',[0.69,0.19,0.38],'linestyle','--','linewidth',2.3);
hold on
plot(Time,varEKFFinal(1,2:end),'color',[0.00,0.79,0.14],'linestyle','--','linewidth',2)
hold on
ylabel('RMSE (m)')
xlabel('Time (s)')
legend('CRB-KF','SCKF','SRCKF','CKF','EKF','CRB','CRB-SCKF','var-EKF')
title('Velocity Estimation')
set(gca,'FontSize',13);
grid on
hold on

subplot('Position',[0.085,0.091644204851752,0.88,0.35])
plot(Time,devLocUpdateFinal(2,2:end),'b','linewidth',2.3);
hold on
plot(Time,devLocSCKF_Final(2,2:end),'color',[0.69,0.19,0.38],'linewidth',2.3);
hold on
plot(Time,devLocSRCKF_Final(2,2:end),'m','linewidth',4);
hold on
plot(Time,devLocCubKF_Final(2,2:end),'g','linewidth',3);
hold on
plot(Time,devLocKal_Final(2,2:end),'k','linestyle','-.','linewidth',2.3);
hold on
plot(Time,real(crbFinal(2,2:end)),'r','linestyle','--','linewidth',2.3);
hold on
plot(Time,real(crbFinalSCKF(2,2:end)),'color',[0.69,0.19,0.38],'linestyle','--','linewidth',2.3);
hold on
plot(Time,varEKFFinal(2,2:end),'color',[0.00,0.79,0.14],'linestyle','--','linewidth',2)
% hold on
ylabel('RMSE (m)')
xlabel('Time (s)')
legend('CRB-KF','SCKF','SRCKF','CKF','EKF','CRB','CRB-SCKF','var-EKF')
title('Position Estimation')
set(gca,'FontSize',13);
grid on
currentDate = date;
name = ['Dat=' currentDate ';a=' num2str(abs(phiMeanInit1(1))) '_' num2str(abs(phiMeanInit1(2))) '_' num2str(abs(phiMeanInit1(3)))...
    ';p=' num2str(Prob1) '_' num2str(Prob2) ...
    ';gridDel=' num2str(delPosAP) ';varMea=' num2str(varDis) ';DeltaT=' num2str(timDelta) ...
    ';MonTimes=' num2str(monTimesReal) '.fig'];
saveas(gcf,name);

% figure
% iniIdx = 1;
% endIdx = length(nodeLocUpdate(dim+1,1:end-1));
% plot3(nodeLocUpdate(dim+1,iniIdx:endIdx),nodeLocUpdate(dim+2,iniIdx:endIdx),nodeLocUpdate(dim+3,iniIdx:endIdx),'c','linewidth',1.5)
% hold on
% plot3(nodeLocEst(1,iniIdx:endIdx),nodeLocEst(2,iniIdx:endIdx),nodeLocEst(3,iniIdx:endIdx),'b','linewidth',1.5)
% hold on
% % plot3(statEstSCKF(dim+1,1:end-1),statEstSCKF(dim+2,1:end-1),statEstSCKF(dim+3,1:end-1),'m','linewidth',1.5)
% % hold on
% plot3(statEstEKF(dim+1,iniIdx:endIdx),statEstEKF(dim+2,iniIdx:endIdx),statEstEKF(dim+3,iniIdx:endIdx),'color',[0.00,0.79,0.14],'linewidth',1.5)
% hold on
% plot3(nodeLocEst1(1,iniIdx:endIdx),nodeLocEst1(2,iniIdx:endIdx),nodeLocEst1(3,iniIdx:endIdx),'r','linewidth',1.5)
% hold on
% % plot3(nodeEstPF(dim+1,1:end-1),nodeEstPF(dim+2,1:end-1),nodeEstPF(dim+3,1:end-1),'c','linewidth',1.5)
% % hold on
% plot3(nodeLocUpdate(dim+1,1),nodeLocUpdate(dim+2,1),nodeLocUpdate(dim+3,1),'bo')
% hold on
% plot3(nodeLocUpdate(dim+1,end-1),nodeLocUpdate(dim+2,end-1),nodeLocUpdate(dim+3,end-1),'r*')
% hold on
% ylabel('dim_2 (m)')
% xlabel('dim_1 (m)')
% legend('Ground Truth','CKF','EKF','Cub-KF','Initial Point','Destination')
% title('Trace of the Target')
% grid on