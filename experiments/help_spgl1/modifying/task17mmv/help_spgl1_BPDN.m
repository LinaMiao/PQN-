cd ../../../../pqnl1
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task17mmv/
clear;

m = 50; n = 128;                      % Measurement matrix is m x n
k = 14;                               % Set sparsity level x0
A = randn(m,n);                       % Random encoding matrix
[Q,R] = qr(A',0);  A = Q';


% group sparse X0
p = 2; nn = n/2;
X0 = zeros(nn,p);
pp = randperm(nn); pp = pp(1:k);
X0(pp,:) = 1e-3*randn(k,p);
B = A * vec(X0) + 0.005 * randn(m,1);
b = B(:);

groups = p;


options.project     = @(x,weight,tau) NormL12_project(groups,x,weight,tau);
options.primal_norm = @(x,weight    ) NormL12_primal(groups,x,weight);
options.dual_norm   = @(x,weight    ) NormL12_dual(groups,x,weight);

% cd ../task7
% addpath(genpath(pwd))
% cd ../task8
% tau = mixNorm(X0,1,2);

tau = 0;
sigma = 1e-3;
[x_spg,r_spg,g_spg,info_spg] = spgl1(A,B(:),tau,sigma,zeros(size(A,2),1),options);
[x_pqn1,r_pqn1,g_pqn1,info_pqn1] = pqnl1_2(A,B(:),tau,sigma,zeros(size(A,2),1),options);

figure;
subplot(3,1,1); plot(vec(X0));title('X0')
subplot(3,1,2); plot(vec(x_spg)); title('x_spg')
subplot(3,1,3); plot(vec(x_pqn1)); title('x_pqn')
title('Multiple Measurement Vector Basis Pursuit');
info_spg
info_pqn1

%% show result
figure('Name','Solution paths')
plot(info_spg.xNorm1,info_spg.rNorm2,info_pqn1.xNorm1,info_pqn1.rNorm2);hold on
scatter(info_spg.xNorm1,info_spg.rNorm2);
scatter(info_pqn1.xNorm1,info_pqn1.rNorm2);hold off
legend('SPGL1_sasha','PQNl1')
axis tight

