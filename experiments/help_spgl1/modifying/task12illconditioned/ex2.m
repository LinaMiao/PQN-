%% addpath for PQN working
%addpath(genpath('/Volumes/Users/linamiao/Dropbox/PQN/'))
cd ../../../../pqnl1;
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task12illconditioned
addpath(genpath(pwd))
rmpath('/Volumes/Users/linamiao/Dropbox/PQN/pqnl1/minConF/')

%stream = RandStream.getGlobalStream;
%reset(stream);


%% sample matrix and options
m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';

opts.decTol = 1e-3;
opts.optTol = 1e-4;
opts.iterations = 100;

% save temp A m n k opts
% clear;
% load temp

%% problem setting
p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
figure;plot(x0)
b  = A*x0;

tau = norm(x0,1);

%% reconstruct
[x_pqn1,r_pqn1,g_pqn1,info_pqn1] = pqnl1_2(A, b, tau, [], zeros(size(x0)), opts);

% using the original minConf_PQN
[x_pqn2,r_pqn2,g_pqn2,info_pqn2] = pqnl1_test(A, b, tau, [], zeros(size(x0)), opts);
