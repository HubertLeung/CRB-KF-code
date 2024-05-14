function du = partial_h(statFore,posAnc,dim)
du = zeros(2,2*dim);
t = 0;
pk = posAnc;
y = statFore(1:dim);
x = statFore(dim+1:2*dim);
s = x - pk;
yn = y;


dhx_1 = s/norm(s);
dhx_2 = yn/norm(s) - s'*yn*s/norm(s)^3;


dhy_1 = t*s/norm(s);
dhy_2 = (s+t*yn)/norm(s) - s'*yn*s*t/norm(s)^3;

du(1,:) = [dhy_1' dhx_1']; 
du(2,:) = [dhy_2' dhx_2']; 




