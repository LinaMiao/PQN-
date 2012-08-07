%% addpath for PQN working
%addpath(genpath('/Volumes/Users/linamiao/Dropbox/PQN/'))
cd ../../../../pqnl1;
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task11lasso

%stream = RandStream.getGlobalStream;
%reset(stream);


%% sample matrix and options
m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';

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
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, b, tau, [], zeros(size(x0)), opts); 
[x_pqn,r_pqn,g_pqn,info_pqn] = pqnl1_2(A, b, tau, [], zeros(size(x0)), opts);

figure;
subplot(2,1,1); plot(x_spg);title('x_spg')
subplot(2,1,2); plot(x_pqn);title('x_pqn')

which spgl1
which pqnl1_2
which minConf_PQN_pqnl1