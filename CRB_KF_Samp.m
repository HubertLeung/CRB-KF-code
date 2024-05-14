function [nodeLoc_pre,crb] = CRB_KF_Samp(measTol,nodeLoc_pre,posAnc,dim,varMeas,nodeLoc_last,R_tol,M_tol,S_tol,statTrans,statIn,accMean,Prob,meaIdx,numSamp)
numIter = 5e2;      % Number of iterations
minNum = 1e-6;      % Newton's iteration threshold
for nn = 1:1:numIter
    [df,Hes] = Hessian(measTol,nodeLoc_pre,posAnc,dim,varMeas,nodeLoc_last,R_tol(:,:,meaIdx),M_tol(:,:,meaIdx),S_tol(:,:,meaIdx),statTrans,statIn,accMean(:,meaIdx),Prob);
    Hes = (Hes+Hes')/2;
    [V,S] = eig(Hes);
    ev = diag(S);
    V1 = V(:,ev>0);
    e1 = ev(ev>0);
    u_nt = -V1*((V1'*df)./e1);

    
    step = 1;
    nodeLoc_pre = nodeLoc_pre + step*u_nt;

    if u_nt'*u_nt < minNum 
        break
    end
    if nn == numIter
        fprintf('CKF does not converge\n');
        monIdx
        meaIdx
    end
end
[crb,~] = CRBSamp(nodeLoc_pre,posAnc,varMeas,nodeLoc_last,statTrans,statIn,accMean(:,meaIdx),R_tol(:,:,meaIdx),M_tol(:,:,meaIdx),S_tol(:,:,meaIdx),Prob,dim,numSamp);