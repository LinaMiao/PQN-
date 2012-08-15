clear;close all
%% addpath for PQN working
cd ../../../../functions;
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task16bpdn/

%stream = RandStream.getGlobalStream;
%reset(stream);

% %problem setting
m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
b  = A*x0;


%% lasso 
tau = norm(x0,1);
opts.optTol = 1e-4;
opts.fid = fopen('spg_lasso.txt','w');
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, b, tau, [], zeros(size(A,2),1), opts); % Find BP sol'n.
opts.fid = fopen('pqn_lasso.txt','w');
opts.optTol = info_spg.rNorm2(end);
[x_pqn1,r_pqn1,g_pqn1,info_pqn1] = pqnl1_2(A, b, tau, [], zeros(size(A,2),1), opts); % Find BP sol'n.
h = figure; 
subplot(2,1,1);plot(x_spg);axis tight;
subplot(2,1,2);plot(x_pqn1);axis tight;
saveas(h,'lasso result.jpg')
info_spg
info_pqn1

save info_lasso info_spg info_pqn1

%% bpdn
b  = A*x0 + 1e-3*rand(size(A,1),1);

opts.fid = fopen('spg_bpdn','w');
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, b, 0, 1e-3, zeros(size(A,2),1), opts); % Find BP sol'n.
sigma_ref = info_spg.rNorm;
opts.fid = fopen('pqn_bpdn','w');
[x_pqn1,r_pqn1,g_pqn1,info_pqn1] = pqnl1_2(A, b, 0, 1e-3, zeros(size(A,2),1), opts,sigma_ref); % Find BP sol'n.
h = figure; 
subplot(2,1,1);plot(x_spg);axis tight;
subplot(2,1,2);plot(x_pqn1);axis tight;
saveas(h,'bpdn result.jpg');
info_spg
info_pqn1

save info_bpdn info_spg info_pqn1

%% show result
h = figure('Name','Solution paths');
plot(info_spg.xNorm1,info_spg.rNorm2,info_pqn1.xNorm1,info_pqn1.rNorm2);hold on
scatter(info_spg.xNorm1,info_spg.rNorm2);
scatter(info_pqn1.xNorm1,info_pqn1.rNorm2);hold off
legend('SPGL1','PQNl1')
axis tight
saveas(h,'solution path.jpg');

