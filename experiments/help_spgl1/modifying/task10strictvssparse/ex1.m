% this experiment is to test whether pqnl1 can work for the expqnle given
% by help spgl1

%% addpath for PQN working
%addpath(genpath('/Volumes/Users/linamiao/Dropbox/PQN/'))
cd ../../../../pqnl1;
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task10strictvssparse

%stream = RandStream.getGlobalStream;
%reset(stream);


%% sample matrix
m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';

opts.decTol = 1e-3;
opts.optTol = 1e-4;
opts.iterations = 100;
opts.nPrevVals = 1; % opt out the nonmonotone line search 
% 
% save temp A m n k opts
% clear;
% load temp

%% problem setting
% strict problem setting
p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
figure;plot(x0)
b0  = A*x0 + 0.005 * randn(m,1);

% compressible problem setting
nn = linspace(0,1,n);
x0_compress = exp(-nn.^.1);
x0_compress = x0_compress - min(x0_compress);
figure;plot(x0_compress)
x0_compress = x0_compress(:);
b_compress  = A*x0_compress + 0.005 * randn(m,1);

%% reconstruct
[x_sparse,r_sparse,g_sparse,info_sparse] = pqnl1_2(A, b0, 0, 1e-3, zeros(size(A,2),1), opts); % Find BP sol'n.
[x_compress,r_compress,g_compress,info_compress] = pqnl1_2(A, b_compress, 0, 1e-3, zeros(size(A,2),1), opts); % Find BP sol'n.
figure('Name','pqn'); 
subplot(2,1,1);plot(x_sparse);subplot(2,1,2);plot(x_compress);


[x_spg1,r_spgl,g_spgl,info_spg1] = spgl1(A, b0, 0, 1e-3, zeros(size(A,2),1), opts);
[x_spg2,r_spg2,g_spg2,info_spg2] = spgl1(A, b_compress, 0, 1e-3, zeros(size(A,2),1), opts); % Find BP sol'n.
figure('Name','spg'); 
subplot(2,1,1);plot(x_spg1);subplot(2,1,2);plot(x_spg2);


%% show result
info_sparse
info_spg1
info_compress
info_spg2

figure('Name','strict sparse Solution paths')
plot(info_sparse.xNorm1,info_sparse.rNorm2,info_spg1.xNorm1,info_spg1.rNorm2);hold on
scatter(info_sparse.xNorm1,info_sparse.rNorm2);
scatter(info_spg1.xNorm1,info_spg1.rNorm2);hold off
legend('pqn','spg')
axis tight

figure('Name','compress signal Solution paths')
plot(info_compress.xNorm1,info_compress.rNorm2,info_spg2.xNorm1,info_spg2.rNorm2);hold on
scatter(info_compress.xNorm1,info_compress.rNorm2);
scatter(info_spg2.xNorm1,info_spg2.rNorm2);hold off
legend('pqn','spg')
axis tight



