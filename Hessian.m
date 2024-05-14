function [df,Hes] = Hessian(measTol,nodeLoc_pre,posAnc,dim,varMeas,nodeLoc_last,R_tol,M_tol,S_tol,statTrans,statIn,accMean,Prob)
dk = measTol(1);
vk = measTol(2);
t = 0;
Dk = varMeas(1);
Dkk = varMeas(2);
Q = varMeas.*eye(2);
pk = posAnc;
y = nodeLoc_pre(1:dim);
x = nodeLoc_pre(dim+1:2*dim);
s = x - pk;
yn = y;

dhx_1 = s/norm(s);
dhx_2 = (yn/norm(s) - s'*yn*s/norm(s)^3);
dhx = [dhx_1 dhx_2];


dhy_1 = t*s/norm(s);
dhy_2 = (s+t*yn)/norm(s) - s'*yn*s*t/norm(s)^3;
dhy = [dhy_1 dhy_2];

dh = [dhy;dhx];


ddfxx_1 = -2/Dk*((dk/norm(s)-1)*eye(dim) - dk*s*(s)'/norm(s)^3);

G1 = -vk*yn*(s)'/norm(s)^3;
G2 = -vk*((s*yn'+(s)'*yn*eye(dim))/(norm(s))^3 - 3*(s)'*yn*s*(s)'/(norm(s))^5);
G3 = -1*(yn*yn')/norm(s)^2 + 2*(s)'*yn*yn*(s)'/norm(s)^4;
G4 = 2*(s)'*yn*s*yn'/norm(s)^4 + ((s)'*yn)^2*eye(dim)/norm(s)^4 - 4*((s)'*yn)^2*(s*(s)')/norm(s)^6;
ddfxx_2 = -2/Dkk*(G1+G2+G3+G4);

ddfxx = ddfxx_1+ddfxx_2;


ddfyy_1 = -2/Dk*((dk/norm(s)-1)*t^2*eye(dim) - dk*t^2*s*(s)'/norm(s)^3);

S1 = 2*vk*t*eye(dim)/norm(s) - vk*t*(s+t*yn)*(s)'/norm(s)^3;
S2 = -1*((s+t*yn)*(s+t*yn)'+ 2*t*(s)'*yn*eye(dim))/norm(s)^2 + 2*t*(s)'*yn*(s+t*yn)*(s)'/norm(s)^4;
S3 = t*(2*(s)'*yn*s*(s+t*yn)'+t*((s)'*yn)^2*eye(dim))/norm(s)^4 - 4*t^2*((s)'*yn)^2*s*(s)'/norm(s)^6;
S4 = -1*vk*t*((s*(s+t*yn)'+t*(s)'*yn*eye(dim))/norm(s)^3 - 3*t*(s)'*yn*(s*s')/norm(s)^5);
ddfyy_2 = -2/Dkk*(S1+S2+S3+S4);

ddfyy = ddfyy_1 + ddfyy_2;


ddfyx_1 = -2/Dk*((dk/norm(s)-1)*t*eye(dim)-dk*t*(s*s')/norm(s)^3);

N1 = vk/norm(s)*(eye(dim) - (s+t*yn)*s'/norm(s)^2);
N2 = -1*vk*t/norm(s)^3*(s*yn'+s'*yn*eye(dim)-3*s'*yn*(s*s')/norm(s)^2);
N3 = -1/norm(s)^2*((s+t*yn)*yn'+s'*yn*eye(dim))+2/norm(s)^4*s'*yn*(s+t*yn)*s';
N4 = t/norm(s)^4*(2*s'*yn*s*yn'+(s'*yn)^2*eye(dim)-4*(s'*yn)^2*(s*s')/norm(s)^2);
ddfyx_2 = -2/Dkk*(N1+N2+N3+N4);

ddfyx = ddfyx_1+ddfyx_2;

ddfxy = ddfyx';

hstep = [norm(s);s'*yn/norm(s)];


df = - 2*dh*(Q\(measTol - hstep));
Hes = [ddfyy ddfyx;ddfxy ddfxx];

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
u1 = -1/2*(nodeLoc_pre - mu1)'*(R_tol\(nodeLoc_pre - mu1));
u2 = -1/2*(nodeLoc_pre - mu2)'*(M_tol\(nodeLoc_pre - mu2));
u3 = -1/2*(nodeLoc_pre - mu3)'*(S_tol\(nodeLoc_pre - mu3));


g1 = -R_tol\(nodeLoc_pre - mu1);
g2 = -M_tol\(nodeLoc_pre - mu2);
g3 = -S_tol\(nodeLoc_pre - mu3);


df_2 = 1/(N1+N2*exp(u2-u1)+N3*exp(u3-u1))*(N1*g1+N2*exp(u2-u1)*g2+N3*exp(u3-u1)*g3);
Hes_2 = -1/(N1+N2*exp(u2-u1)+N3*exp(u3-u1))^2*(N1*g1 + N2*exp(u2-u1)*g2 + N3*exp(u3-u1)*g3)*(N1*g1 + N2*exp(u2-u1)*g2 + N3*exp(u3-u1)*g3)' + ...
    1/(N1+N2*exp(u2-u1)+N3*exp(u3-u1)) * (N1*(g1*g1') - N1*(R_tol\eye(numElement)) + N2*exp(u2-u1)*(g2*g2') - N2*exp(u2-u1)*(M_tol\eye(numElement))...
    + N3*exp(u3-u1)*(g3*g3') - N3*exp(u3-u1)*(S_tol\eye(numElement)));


df = real(df - df_2);
Hes = real(Hes - Hes_2);





