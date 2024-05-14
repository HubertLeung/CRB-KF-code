function [crb,flag] = CRBSamp(nodeLoc_pre,posAnc,varMeas,nodeLoc_last,statTrans,statIn,accMean,R_tol,M_tol,S_tol,Prob,dim,numSamp)
numElement = 2*dim;
du = zeros(2,numElement);
Q = diag(varMeas);
t = 0;
pk = posAnc;
y = nodeLoc_pre(1:dim);
x = nodeLoc_pre(dim+1:2*dim);
s = x - pk; 

dhx_1 = s/norm(s);
dhx_2 = y/norm(s) - s'*y*s/norm(s)^3;

dhy_1 = t*s/norm(s);
dhy_2 = (s+t*y)/norm(s) - s'*y*s*t/norm(s)^3;

du(1,:) = -[dhy_1' dhx_1']; 
du(2,:) = -[dhy_2' dhx_2']; 

F1 = du'*(Q\du);
%%
q = zeros(numElement,numElement);
Prob1 = Prob(1);
Prob2 = Prob(2);
Prob3 = Prob(3);

accMean1 = accMean(1:dim);
accMean2 = accMean(dim+1:2*dim);
accMean3 = accMean(2*dim+1:3*dim);

mu1 = statTrans*nodeLoc_last + statIn*accMean1;
mu2 = statTrans*nodeLoc_last + statIn*accMean2;
mu3 = statTrans*nodeLoc_last + statIn*accMean3;

N1 = Prob1/sqrt((2*pi)^numElement*det(R_tol));
N2 = Prob2/sqrt((2*pi)^numElement*det(M_tol));
N3 = Prob3/sqrt((2*pi)^numElement*det(S_tol));

numReal = numSamp;
data = zeros(numElement,numSamp);
noise = zeros(numElement,numSamp);
for j = 2:1:numSamp+1
    chi = randsrc(numElement,1,[1,2,3;Prob1,Prob2,Prob3]);
    pdf1 = statIn*accMean1 + real(sqrtm(R_tol))*randn(numElement,1);
    pdf2 = statIn*accMean2 + real(sqrtm(M_tol))*randn(numElement,1);
    pdf3 = statIn*accMean3 + real(sqrtm(S_tol))*randn(numElement,1);
    noise(:,j-1) = pdf1.*(chi==1)+pdf2.*(chi==2)+pdf3.*(chi==3);
    data(:,j-1) = statTrans*nodeLoc_last + noise(:,j-1);
end

for i = 2:1:numSamp+1
    nodeLoc = data(:,i-1); 
    f1 = N1*exp(-1/2*(nodeLoc - mu1)'*(R_tol\(nodeLoc - mu1)));
    f2 = N2*exp(-1/2*(nodeLoc - mu2)'*(M_tol\(nodeLoc - mu2)));
    f3 = N3*exp(-1/2*(nodeLoc - mu3)'*(S_tol\(nodeLoc - mu3)));

    g1 = -R_tol\(nodeLoc - mu1);
    g2 = -M_tol\(nodeLoc - mu2);
    g3 = -S_tol\(nodeLoc - mu3);

    u1 = -1/2*(nodeLoc - mu1)'*(R_tol\(nodeLoc - mu1));
    u2 = -1/2*(nodeLoc - mu2)'*(M_tol\(nodeLoc - mu2));
    u3 = -1/2*(nodeLoc - mu3)'*(S_tol\(nodeLoc - mu3));


    q_step_1 = -1/(N1+N2*exp(u2-u1)+N3*exp(u3-u1))^2*(N1*g1 + N2*exp(u2-u1)*g2 + N3*exp(u3-u1)*g3)*(N1*g1 + N2*exp(u2-u1)*g2 + N3*exp(u3-u1)*g3)' + ...
        1/(N1+N2*exp(u2-u1)+N3*exp(u3-u1)) * (N1*(g1*g1') - N1*(R_tol\eye(numElement)) + N2*exp(u2-u1)*(g2*g2') - N2*exp(u2-u1)*(M_tol\eye(numElement))...
        + N3*exp(u3-u1)*(g3*g3') - N3*exp(u3-u1)*(S_tol\eye(numElement)));
    
    q = q + q_step_1;
end
F2 = -q/numReal;
F = F1 + F2;
flag = 0;
if sum(isnan(F),'all') > 0 || sum(isinf(F),'all') > 0
    flag = 1;
end

if isnan(F) 
    1;
    crb = 2*eye(numElement);
else
    crb = real(F\eye(numElement));
end

