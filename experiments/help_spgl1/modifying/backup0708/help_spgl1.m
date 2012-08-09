% this experiment is to test whether pqnl1 can work for the expqnle given
% by help spgl1

%% addpath for PQN working
%addpath(genpath('/Volumes/Users/linamiao/Dropbox/PQN/'))
cd ../../../..;
addpath(genpath(pwd))
cd ./experiments/help_spgl1/modifying/backup0708

%stream = RandStream.getGlobalStream;
%reset(stream);

%problem setting
m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
b  = A*x0 + 0.005 * randn(m,1);

opts.decTol = 1e-3;
opts.optTol = 1e-4;
opts.iterations = 100;
opts.nPrevVals = 1; % opt out the nonmonotone line search 

% save temp A b opts
% clear
% load temp.mat
%% spgl1
tic
%[x_spg,r_spg,g_spg,info_spg] = spgl1(A, b, 0, 1e-3, [], opts); % Find BP sol'n.
[x_spg,r_spg,g_spg,info_spg] = spgl1_origin(A, b, 0, 1e-3, [], opts); % Find BP sol'n.
toc
%% pqnl1
%opts.iterations = .25*opts.iterations;
[x_pqn,r_pqn,g_pqn,info_pqn] = pqnl1_2(A, b, 0, 1e-3, zeros(size(A,2),1), opts); % Find BP sol'n.
toc
%% show result
figure; subplot(2,1,1);plot(x_spg);subplot(2,1,2);plot(x_pqn);
info_spg
info_pqn

%% draw solution path
figure('Name','Solution paths')
plot(info_spg.xNorm1,info_spg.rNorm2,info_pqn.xNorm1,info_pqn.rNorm2);hold on
scatter(info_spg.xNorm1,info_spg.rNorm2);
scatter(info_pqn.xNorm1,info_pqn.rNorm2);hold off
xlabel('one-norm model')
ylabel('two-norm residual')
title('Solutions paths')
legend('SPGl1','PQNl1')
axis tight

%% check functions
% open ./minConF_PQN_2.m
% open ./pqnl1_2.m
