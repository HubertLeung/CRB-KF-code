function [stateUpd, PSR, root_PSR] = SRCKF(varMeas,dim,posAnc,measTol,delAccMean,delAccVar,statTrans,statIn,nodeLoc_last1,root_PSR,numElement)
varphi = sqrt(numElement)*[eye(numElement) -eye(numElement)];
% Prediction

cubaturePointNum = 2*numElement;
choleskyDecompose = root_PSR;

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

delta_Pre = 1/sqrt(cubaturePointNum)*(statePre_CP - kron(statePre,ones(1,length(statePre_CP))));
Sq = statIn*chol(delAccVar,'lower');

[~,rootCovarPre] = qr([delta_Pre Sq]',0);
rootCovarPre = rootCovarPre';

% Update
stateUpdate_CP = zeros(numElement,cubaturePointNum);
for iCubaturePoint = 1:cubaturePointNum
    stateUpdate_CP(:,iCubaturePoint) = statePre + rootCovarPre*varphi(:,iCubaturePoint);
end

measPre_CP = zeros(length(measTol),cubaturePointNum);
for idx = 1:cubaturePointNum
    measPre_CP(:,idx) = mea_trans(stateUpdate_CP(:,idx),posAnc,dim);
end
measPre = measPre_CP*weight_CP';

delta_measUp = 1/sqrt(cubaturePointNum)*(measPre_CP - kron(measPre,ones(1,length(statePre_CP))));
Sr = sqrtm(diag(varMeas));
[~,rootMeaUpzz] = qr([delta_measUp Sr]',0);
rootMeaUpzz = rootMeaUpzz';

delta_covUp = 1/sqrt(cubaturePointNum)*(stateUpdate_CP - kron(statePre,ones(1,length(statePre_CP))));

rootCovUpxz = delta_covUp*delta_measUp';

Krs = (rootCovUpxz/rootMeaUpzz')/rootMeaUpzz;
stateUpd = statePre + Krs*(measTol-measPre);

[~,root_PSR] = qr([delta_covUp-Krs*delta_measUp Krs*Sr]',0);
root_PSR = root_PSR';
PSR = root_PSR*root_PSR';