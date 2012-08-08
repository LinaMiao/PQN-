% this experiment is to test whether pqnl1 can work for the expqnle given
% by help spgl1

%% addpath for PQN working
%addpath(genpath('/Volumes/Users/linamiao/Dropbox/PQN/'))
cd ../../../../pqnl1;
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task16bpdn

%stream = RandStream.getGlobalStream;
%reset(stream);

% %problem setting
m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
b  = A*x0;
% 
% % opts.decTol = 1e-3;
% % opts.optTol = 1e-4; 
% opts.iterations = 500;
% opts.nPrevVals = 1; % opt out the nonmonotone line search 
% 
% save temp A b x0 opts
clear
load temp.mat


%% lasso 
tau = norm(x0,1);
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, b, tau, [], zeros(size(A,2),1), opts); % Find BP sol'n.

[x_pqn1,r_pqn1,g_pqn1,info_pqn1] = pqnl1_2(A, b, tau, [], zeros(size(A,2),1), opts); % Find BP sol'n.

figure; subplot(2,1,1);plot(x_spg);subplot(2,1,2);plot(x_pqn1);
info_spg
info_pqn1

%% show result
figure('Name','Solution paths')
plot(info_spg.xNorm1,info_spg.rNorm2,info_pqn1.xNorm1,info_pqn1.rNorm2);hold on
scatter(info_spg.xNorm1,info_spg.rNorm2);
scatter(info_pqn1.xNorm1,info_pqn1.rNorm2);hold off
legend('SPGL1_sasha','PQNl1')
axis tight

%% bpdn
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, b, 0, 1e-3, zeros(size(A,2),1), opts); % Find BP sol'n.

[x_pqn1,r_pqn1,g_pqn1,info_pqn1] = pqnl1_2(A, b, 0, 1e-3, zeros(size(A,2),1), opts); % Find BP sol'n.

figure; subplot(2,1,1);plot(x_spg);subplot(2,1,2);plot(x_pqn1);
info_spg
info_pqn1

%% show result
figure('Name','Solution paths')
plot(info_spg.xNorm1,info_spg.rNorm2,info_pqn1.xNorm1,info_pqn1.rNorm2);hold on
scatter(info_spg.xNorm1,info_spg.rNorm2);
scatter(info_pqn1.xNorm1,info_pqn1.rNorm2);hold off
legend('SPGL1_sasha','PQNl1')
axis tight

%% check functions
% open ./minConF_PQN_2.m
% open ./pqnl1_2.m
