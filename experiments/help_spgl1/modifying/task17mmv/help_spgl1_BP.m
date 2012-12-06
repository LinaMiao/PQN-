clear;close all;

cd ../../../../functions/
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task17mmv/

m = 120; n = 512;                      % Measurement matrix is m x n
k = 10;                               % Set sparsity level x0
A = randn(m,n);                       % Random encoding matrix
[Q,R] = qr(A',0);  A = Q';


% group sparse X0
p = 4; nn = n/p;
X0 = zeros(nn,p);
pp = randperm(nn); pp = pp(1:k);
X0(pp,:) = sign(randn(k,p));
B = A * vec(X0);% + 1e-3*rand(size(A,1),1);
b = B(:);

groups = p;

%save temp A B p X0 groups
load temp
options.project     = @(x,weight,tau) NormL12_project(groups,x,weight,tau);
options.primal_norm = @(x,weight    ) NormL12_primal(groups,x,weight);
options.dual_norm   = @(x,weight    ) NormL12_dual(groups,x,weight);


tau = 0;
sigma = 1e-3;
options.iterations = 800;
options.fid = fopen('spg','w');
[x_spg,r_spg,g_spg,info_spg] = spgl1(A,B(:),0,sigma,zeros(size(A,2),1),options);
options.fid = fopen('pqn','w');
sigma_ref = info_spg.rNorm;
[x_pqn1,r_pqn1,g_pqn1,info_pqn1] = pqnl1_2(A,B(:),0,sigma,zeros(size(A,2),1),options,sigma_ref);


h = figure;
subplot(3,1,1); plot(vec(X0));title('X0'); ylim([-1 1])
subplot(3,1,2); plot(vec(x_spg)); title('spg recovery, 228+170 mat-vec'); ylim([-1 1])
subplot(3,1,3); plot(vec(x_pqn1)); title('pqn recovery, 100+100 mat-vec'); ylim([-1 1])
saveas(h,'mmv');

info_spg
info_pqn1

save info info_spg info_pqn1
%% show result
% show gourpsparsity of X0
h = figure;
spy(X0);pbaspect([1 1 1]);
saveas(h,'spyofx0')
% solution path
h = figure('Name','Solution paths');
plot(info_spg.xNorm1,info_spg.rNorm2,info_pqn1.xNorm1,info_pqn1.rNorm2);hold on
scatter(info_spg.xNorm1,info_spg.rNorm2);
scatter(info_pqn1.xNorm1,info_pqn1.rNorm2);hold off
legend('SPGL1_sasha','PQNl1')
axis tight
saveas(h,'solution path');
