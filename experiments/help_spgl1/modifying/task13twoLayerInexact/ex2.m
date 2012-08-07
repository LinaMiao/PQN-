%% addpath for PQN working
%addpath(genpath('/Volumes/Users/linamiao/Dropbox/PQN/'))
cd ../../../../pqnl1;
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task13twoLayerInexact
addpath(genpath(pwd))
rmpath('/Volumes/Users/linamiao/Dropbox/PQN/pqnl1/minConF/')

%stream = RandStream.getGlobalStream;
%reset(stream);


%% sample matrix and options
% m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
% A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
% 
% opts.decTol = 1e-3;
% opts.optTol = 1e-4;
% opts.iterations = 100;
% 
% p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
% figure;plot(x0)
% b  = A*x0;
% 
% tau = norm(x0,1);
% 
% save temp A m n k b tau x0 opts
clear;close all;
load temp


%% reconstruct
%[x_spg1,r_spg1,g_spg1,info_spg1] = spgl1(A, b, tau, [], zeros(size(x0)), opts);
flag = 1;
[x_pqn1,r_pqn1,g_pqn1,info_pqn1] = pqnl1_2(A, b, tau, [], zeros(size(x0)), opts,flag);
flag = 0;
[x_pqn2,r_pqn2,g_pqn2,info_pqn2] = pqnl1_2(A, b, tau, [], zeros(size(x0)), opts,flag);



figure;
subplot(3,1,1);plot(x0);title('x0')
subplot(3,1,2);plot(x_pqn1);title('two layer')
subplot(3,1,3);plot(x_pqn2);title('optTol')

%% spgl1 as a reference
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, b, tau, [], zeros(size(x0)), opts);

info_pqn1
info_pqn2
info_spg

