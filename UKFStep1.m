function [xhatFinal, PFinal] = UKFStep1(xhat, zk, P, Q, R, a, k, b,dim,Prob,accMean,statTrans,statIn,posAnc,numVair,numPoc)
L = 2*dim;
% Init - Create the weights
[Wa, Wc] = ukf_create_weights(a, b, k, L);

% PREDICT: Step 0 - Predict the state and state estimation error covariance at the next time step
s = ukf_compute_sigma_points(xhat, P, a, k, L);

% PREDICT: Step 1 - Run our transition function
x = ukf_transition(s, L, Prob,accMean,statTrans,statIn,dim);
% PREDICT: Step 2 - Combine the predicted states to obtain the predicted states
xhat = ukf_multiply_weights(x, Wa, L);

% PREDICT: Step 3 - Compute the covariance of the predicted state
P = ukf_estimate_covariance(x, xhat, Wc, Q, L);

h = partial_h(xhat,posAnc,dim);
K = P*h'/(h*P*h'+ kron(R,eye(numPoc)));
xhatFinal = xhat + K*(zk - mea_trans(xhat,posAnc,dim));
PFinal = (eye(2*dim)-K*h)*P*(eye(2*dim)-K*h)'+ K*R*K';

end

function [Wa, Wc] = ukf_create_weights(a, b, k, L)
N = 2 * L + 1; 
Wa = zeros(1, N); 
Wc = zeros(1, N); 
for i = 1:N
    if(i == 1) 
        Wa(i) = (a*a*k-L)/(a*a*k); 
        Wc(i) = Wa(i) + 1 - a*a + b;
    else
        Wa(i) = 1/(2*a*a*k);
        Wc(i) = Wa(i);
    end
end
end

function [s] = ukf_compute_sigma_points(x, P, a, k, L)
N = 2 * L + 1;
compensate = L + 1;
s = zeros(L, N);
1;
[V,D] = eig(P);
[idx1,idx2] = find(D<0);

if isempty(idx2)
    1;
else
    D(idx1,idx2) = 0;
end
P = real(V*D*V') + 1e-6*eye(size(P));

A = a*sqrt(k)*chol(P);
for j = 1:N
    if(j == 1)
        s(:, j) = x;
    elseif(and(j >= 2, j <= L + 1))
        s(:, j) = x + A(:, j - 1);
    else
        s(:, j) = x - A(:, j - compensate);
    end
end
end

function x = ukf_transition(s, L, Prob,accMean,statTrans,statIn,dim)
Prob1 = Prob(1);
Prob2 = Prob(2);
Prob3 = Prob(3);
N = 2 * L + 1;
x = zeros(L, N);
for i = 1:N
    x(:,i) = statTrans*s(:,i) + Prob1*statIn*accMean(1:dim) + Prob2*statIn*accMean(dim+1:2*dim) + Prob3*statIn*accMean(2*dim+1:3*dim);
end
end

function x = ukf_multiply_weights(xi, W, L)
N = 2 * L + 1;
x = zeros(L, 1);
for i = 1:N
    x = x + W(i)*xi(:, i);
end
end

function z = ukf_multiply_weights_mea(zi, W, L,numVair)
N = 2 * L + 1;
z = zeros(numVair, 1);
for i = 1:N
    z = z + W(i)*zi(:, i);
end
end

function P = ukf_estimate_covariance(xi, x, W, O, L)
N = 2 * L + 1;
P = zeros(L, L);
for i = 1:N
    P = P + W(i)*(xi(:, i) - x)*(xi(:, i) - x)';
end
P = P + O;
end

function P = ukf_estimate_covariance_mea(zi, z, W, O, L, numVair)
N = 2 * L + 1;
P = zeros(numVair, numVair);
for i = 1:N
    P = P + W(i)*(zi(:, i) - z)*(zi(:, i) - z)';
end
P = P + O;
end

function Csz = ukf_estimate_cross_covariance(s, xhat, z, zhat, Wc, L, numVair)
N = 2 * L + 1;
Csz = zeros(L, numVair);
for i = 1:N
    Csz = Csz + Wc(i)*(s(:, i) - xhat)*(z(:, i) - zhat)';
end
end

function K = ukf_create_kalman_K(Shat, Csz, L, numVair)
K = Csz/Shat;
1;
% K = zeros(numVair,L);
% 1;
% Csz = Csz';
% for i = 1:L
%     % Solve Ax = b with Cholesky
%     A = chol(Shat, 'lower');
%     y = linsolve(A, Csz(:, i));
%     K(:,i) = linsolve(A', y);
% end
% K = K';

end

function [xhat, P] = ukf_state_update(K, Shat, P, xhat, zk, zhat, L)
xhat = xhat + K*(zk - zhat);
P = P - K*Shat*K';
end
