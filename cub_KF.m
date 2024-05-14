function [stateUpd, P1] = cub_KF(varMeas,dim,posAnc,measTol,delAccMean,delAccVar,statTrans,statIn,nodeLoc_last1,P1,numElement)
varphi = sqrt(numElement)*[eye(numElement) -eye(numElement)];
% Prediction
% Generate cubature points
cubaturePointNum = 2*numElement;
choleskyDecompose = chol(P1,'lower');% square root of the covariance matrix

state_CP = zeros(numElement,cubaturePointNum);
for iCubaturePoint = 1:cubaturePointNum
    state_CP(:,iCubaturePoint) = nodeLoc_last1 + choleskyDecompose*varphi(:,iCubaturePoint);
end

weight_CP = 1/cubaturePointNum.*ones(1,cubaturePointNum);


cubaturePoint_num = size(state_CP,2);
statePre_CP = zeros(numElement,cubaturePoint_num);
for idx = 1:cubaturePoint_num
    statePre_CP(:,idx) = statTrans*state_CP(:,idx) + statIn*delAccMean;
end

statePre = statePre_CP*weight_CP';

covarPre = zeros(numElement,numElement);
for idx = 1:cubaturePoint_num
    covarPre = covarPre + statePre_CP(:,idx)*statePre_CP(:,idx)';
end
covarPre = covarPre/cubaturePoint_num - statePre*statePre';

covarPre = covarPre + statIn*delAccVar*statIn';

% Update

choleskyDecompose = chol(covarPre,'lower');

stateUpdate_CP = zeros(numElement,cubaturePointNum);% cubature点状态
for iCubaturePoint = 1:cubaturePointNum
    stateUpdate_CP(:,iCubaturePoint) = statePre + choleskyDecompose*varphi(:,iCubaturePoint);
end

measPre_CP = zeros(length(measTol),cubaturePointNum);
for idx = 1:cubaturePointNum
    measPre_CP(:,idx) = mea_trans(stateUpdate_CP(:,idx),posAnc,dim);
end
measPre = measPre_CP*weight_CP';

Covar = zeros(numElement,length(measTol));
Scov = zeros(length(measTol),length(measTol));
sigmaPoint_num = size(state_CP,2);
for idx = 1:sigmaPoint_num
    Covar = Covar + stateUpdate_CP(:,idx)*measPre_CP(:,idx)';
    Scov = Scov + measPre_CP(:,idx)*measPre_CP(:,idx)';
end
Covar = Covar/sigmaPoint_num - statePre*measPre';
Scov = Scov/sigmaPoint_num - measPre*measPre';

Scov = Scov + diag(varMeas);

K1 = Covar/Scov;
stateUpd = statePre + K1*(measTol-measPre);

P1 = covarPre - K1*Scov*K1';
