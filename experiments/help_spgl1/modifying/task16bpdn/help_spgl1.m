clear;close all
%% addpath for PQN working
cd ../../../../functions;
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task16bpdn/

%stream = RandStream.getGlobalStream;
%reset(stream);
load v
rng(v);
% %problem setting
m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
%p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
%A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
%b  = A*x0;
load temp

%% lasso 
tau = norm(x0,1);
opts.optTol = 1e-4;
opts.fid = fopen('spg_lasso.txt','w');
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, b, tau, [], zeros(size(A,2),1), opts); % Find BP sol'n.
opts.fid = fopen('pqn_lasso.txt','w');
sigma_ref = info_spg.rNorm;
corrections_list = [1 2 3 4 5 10 15 20 25 30];
for i = 1:length(corrections_list)
    corrections = corrections_list(i);
    [x_pqn1(:,i),r_pqn1,g_pqn1,info_pqn1(:,i)] = pqnl1_2(A, b, tau, [], zeros(size(A,2),1), opts, sigma_ref,corrections); % Find BP sol'n.
    h = figure(i); 
    subplot(3,1,1);plot(x0);axis tight;title('true signal')
    subplot(3,1,2);plot(x_spg);axis tight; title('spg recovery')%, 86+67 mat-vec')
    subplot(3,1,3);plot(x_pqn1(:,i));axis tight; title('pqn recovery')%, 38+38 mat-vec')
    saveas(h,'lasso result')
end
info_spg
info_pqn1

h = figure; plot(corrections_list,2.*[info_pqn1.nProdA],'o-')
xlabel('the number of gradient vectors needed to save in approximate Hessian')
ylabel('the number of mat-vec needed')
title('comparision of mat-vec between SPGl1 and PQNl1 in one simple lasso problem')
saveas(h,'comparision of mat-vec between SPGl1 and PQNl1 in one simple lasso problem')
save info_lasso info_spg info_pqn1

% solution path for lasso problem, proving pqn does converge fast, i.e.
% within fewer iterations
h = figure('Name','Solution paths lasso');
plot(info_spg.xNorm1,info_spg.rNorm2,[0;info_pqn1.xNorm1],[info_spg.rNorm2(1);info_pqn1.rNorm2]);hold on
scatter(info_spg.xNorm1,info_spg.rNorm2);
scatter(info_pqn1.xNorm1,info_pqn1.rNorm2);hold off
legend('SPGL1','PQNl1')
axis tight
xlim([0 1.5*info_spg.xNorm1(end)])
saveas(h,'solution path lasso');

%% bpdn
%b  = A*x0 + 1e-3*rand(size(A,1),1);
%save temp1 A b x0
load temp1 
opts.fid = fopen('spg_bpdn.txt','w');
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, b, 0, 1e-3, zeros(size(A,2),1), opts); % Find BP sol'n.
sigma_ref = info_spg.rNorm;
opts.fid = fopen('pqn_bpdn.txt','w');
[x_pqn1,r_pqn1,g_pqn1,info_pqn1] = pqnl1_2(A, b, 0, 1e-3, zeros(size(A,2),1), opts,sigma_ref); % Find BP sol'n.
h = figure; 
subplot(2,1,1);plot(x_spg);axis tight;
subplot(2,1,2);plot(x_pqn1);axis tight;
saveas(h,'bpdn result');
info_spg
info_pqn1

save info_bpdn info_spg info_pqn1

%% show result
h = figure('Name','Solution paths bpdn');
plot(info_spg.xNorm1,info_spg.rNorm2,[0;info_pqn1.xNorm1],[info_spg.rNorm2(1);info_pqn1.rNorm2]);hold on
scatter(info_spg.xNorm1,info_spg.rNorm2);
scatter(info_pqn1.xNorm1,info_pqn1.rNorm2);hold off
legend('SPGL1','PQNl1')
axis tight
saveas(h,'solution path bpdn');

