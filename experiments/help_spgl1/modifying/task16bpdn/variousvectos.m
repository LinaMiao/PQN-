% this script is to see what the influence of storing vectors increase
clear;close all
%% addpath for PQN working
cd ../../../../functions;
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task16bpdn/

%stream = RandStream.getGlobalStream;
%reset(stream);
% load v
% rng(v);


totalExperiment = 20;
for nexperiment = 1:totalExperiment

    % %problem setting
    m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
    p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
    A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
    b  = A*x0;


    %% lasso 
    tau = norm(x0,1);
    opts.optTol = 1e-4;
    opts.fid = fopen('spg_lasso.txt','w');
    [x_spg,r_spg,g_spg,info_spg(nexperiment)] = spgl1(A, b, tau, [], zeros(size(A,2),1), opts); % Find BP sol'n.
    opts.fid = fopen('pqn_lasso.txt','w');
    sigma_ref = info_spg(nexperiment).rNorm;
    corrections_list = [1 2 3 4 5 10 15 20 25 30];
    for i = 1:length(corrections_list)
        corrections = corrections_list(i);
        [~,r_pqn1,g_pqn1,info_pqn1(i,nexperiment)] = pqnl1_2(A, b, tau, [], zeros(size(A,2),1), opts, sigma_ref,corrections); % Find BP sol'n.
    end
    save info info_pqn1
end
tempspg = [info_spg.nProdA] + [info_spg.nProdAt];
tempspg = reshape(tempspg,size(info_spg));
tempspg = mean(tempspg) * ones(i,1);

temp = [info_pqn1.nProdA] + [info_pqn1.nProdAt];
temp = reshape(temp,size(info_pqn1));
PQNnProdA = mean(temp,2);

h = figure;
%plot(corrections_list,tempspg,'r');hold on;
plot(corrections_list,PQNnProdA,'-o');
xlabel('the number of gradient vectors needed to save in approximate Hessian')
ylabel('the number of mat-vec needed');
axis tight
gtext({'mat-vec of SPG is 319'})%,num2str(round(tempspg(1)))})
%ylim([50 100])
title('mat-vec of PQNl1 with various rank of LBFGS updata')
saveas(h,'comp_spg_pqn_20')
