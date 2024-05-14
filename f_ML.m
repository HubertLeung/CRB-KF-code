function f = f_ML(measTol,nodeLoc_pre,posAnc,dim,varMeas,nodeLoc_last,R_tol,M_tol,S_tol,statTrans,statIn,accMean,Prob)
Q = varMeas.*eye(2);
pk = posAnc;
y = nodeLoc_pre(1:dim);
x = nodeLoc_pre(dim+1:2*dim);
s = x - pk;
yn = y;
hstep = [norm(s);s'*yn/norm(s)];
% f = -1/2*(measTol - hstep)'*(Q\(measTol - hstep));
f = (measTol - hstep)'*(Q\(measTol - hstep));
%%
numElement = 2*dim;
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

f1 = N1*exp(-1/2*(nodeLoc_pre - mu1)'*(R_tol\(nodeLoc_pre - mu1)));
f2 = N2*exp(-1/2*(nodeLoc_pre - mu2)'*(M_tol\(nodeLoc_pre - mu2)));
f3 = N3*exp(-1/2*(nodeLoc_pre - mu3)'*(S_tol\(nodeLoc_pre - mu3)));

g = -log(1e-20 + f1+f2+f3);
f = real(f + g);


