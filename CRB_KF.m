function [nodeLoc_pre,crb] = CRB_KF(measTol,nodeLoc_pre,posAnc,dim,varMeas,nodeLoc_last,R_tol,M_tol,S_tol,statTrans,statIn,accMean,Prob,meaIdx)
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

    % alpha∈[0.01,0.3],beta∈[0.1,0.8]
    step = 1;
    nodeLoc_delta = nodeLoc_pre + step*u_nt;
%     alpha = 0.03;
%     beta = 0.8;
%     while 1
%         nodeLoc_delta = nodeLoc_pre + step*u_nt;
%         f = f_ML(measTol,nodeLoc_pre,posAnc,dim,varMeas,nodeLoc_last,R_tol(:,:,meaIdx),M_tol(:,:,meaIdx),S_tol(:,:,meaIdx),statTrans,statIn,accMean(:,meaIdx),Prob);
%         f_delta = f_ML(measTol,nodeLoc_delta,posAnc,dim,varMeas,nodeLoc_last,R_tol(:,:,meaIdx),M_tol(:,:,meaIdx),S_tol(:,:,meaIdx),statTrans,statIn,accMean(:,meaIdx),Prob);
%         stt = alpha*step*df'*u_nt;
%         sss = f + stt;
%         if f_delta > sss
%             step = beta*step;
%         else
%             break
%         end
%     end
    nodeLoc_pre = nodeLoc_delta;

    if u_nt'*u_nt < minNum 
        break
    end
    if nn == numIter
        fprintf('CKF does not converge\n');
        monIdx
        meaIdx
    end
end
[crb,~] = CRB(nodeLoc_pre,posAnc,varMeas,nodeLoc_last,statTrans,statIn,accMean(:,meaIdx),R_tol(:,:,meaIdx),M_tol(:,:,meaIdx),S_tol(:,:,meaIdx),Prob,dim);