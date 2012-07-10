% %% generate data if not exist
% if not(exist('camenbert.mat'))
%      MGNFWI_camenbert;
% end

%%
addpath(genpath('../../2DModGaussNewton'))
addpath(genpath('../../..'))

%% inversion

clear
MGNFWI_camenbert;
%load camenbert
opts.iterations = 300;
[results_1] = MGNFWI(reshape(m0,n),Dobs,wavelet,model,opts,1);
save results_1 results_1

clear
MGNFWI_camenbert;
%load camenbert
opts.iterations = 300; 
[results_2] = MGNFWI(reshape(m0,n),Dobs,wavelet,model,opts,2);
save results_2 results_2 % pqnl1

clear
MGNFWI_camenbert;
%load camenbert
opts.iterations = 300;
[results_3] = MGNFWI(reshape(m0,n),Dobs,wavelet,model,opts,3);
save results_3 results_3 % pqnl1_2


%%
open MGNFWI_camenbert
open GN_up
open pqnl1_2
open minconf_pqn_2
open lbfgsupdate