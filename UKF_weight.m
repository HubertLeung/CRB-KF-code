function [stateUpdUKF, PUKF] = UKF_weight(varMeas,dim,posAnc,measTol,delAccMean,delAccVar,statTrans,statIn,nodeLoc_last1,PUKF,numElement)
W0 = 0.09;

varphi = sqrt(numElement/(1-W0))*[eye(numElement) -eye(numElement)];

cubaturePointNum = 2*numElement;
choleskyDecompose = chol(PUKF,'lower');

state_CP_0 = nodeLoc_last1;

state_CP = zeros(numElement,cubaturePointNum); 
for iCubaturePoint = 1:cubaturePointNum
    state_CP(:,iCubaturePoint) = nodeLoc_last1 + choleskyDecompose*varphi(:,iCubaturePoint);
end

weight_CP = (1-W0)/(cubaturePointNum).*ones(1,cubaturePointNum);


cubaturePoint_num = size(state_CP,2);
statePre_CP = zeros(numElement,cubaturePoint_num);
for idx = 1:cubaturePoint_num
    statePre_CP(:,idx) = statTrans*state_CP(:,idx) + statIn*delAccMean;
end

statePre = statePre_CP*weight_CP' + W0*state_CP_0;

covarPre = zeros(numElement,numElement);
for idx = 1:cubaturePoint_num
    covarPre = covarPre + weight_CP(1,idx)*(statePre_CP(:,idx)-statePre)*(statePre_CP(:,idx)-statePre)';
end

covarPre = covarPre + W0*(state_CP_0-statePre)*(state_CP_0-statePre)' + statIn*delAccVar*statIn';


choleskyDecompose = chol(covarPre,'lower');
stateUpdate_CP_0 = statePre;
stateUpdate_CP = zeros(numElement,cubaturePointNum);% cubature点状态
for iCubaturePoint = 1:cubaturePointNum
    stateUpdate_CP(:,iCubaturePoint) = statePre + choleskyDecompose*varphi(:,iCubaturePoint);
end

measPre_CP_0 = mea_trans(stateUpdate_CP_0,posAnc,dim);

measPre_CP = zeros(length(measTol),cubaturePointNum);
for idx = 1:cubaturePointNum
    measPre_CP(:,idx) = mea_trans(stateUpdate_CP(:,idx),posAnc,dim);
end
measPre = measPre_CP*weight_CP' + W0*measPre_CP_0;

Covar = zeros(numElement,length(measTol));
Scov = zeros(length(measTol),length(measTol));
sigmaPoint_num = size(state_CP,2);
for idx = 1:sigmaPoint_num
    Covar = Covar + weight_CP(idx)*(stateUpdate_CP(:,idx)-statePre)*(measPre_CP(:,idx)-measPre)';
    Scov = Scov + weight_CP(idx)*(measPre_CP(:,idx)-measPre)*(measPre_CP(:,idx)-measPre)';
end
Covar = Covar + W0*(stateUpdate_CP_0-statePre)*(measPre_CP_0-measPre)';
Scov = Scov + W0*(measPre_CP_0-measPre)*(measPre_CP_0-measPre)' + diag(varMeas);

KUKF = Covar/Scov;
stateUpdUKF = statePre + KUKF*(measTol-measPre);

PUKF = covarPre - KUKF*Covar';