%%
addpath(genpath('../..'))
MGNFWI_camenbert

%% inversion
clear
load camenbert
opts.iterations = 20;
[results_1] = MGNFWI(reshape(m0,n),Dobs,wavelet,model,opts,1);
save results_1 results_1
clear
load camenbert
[results_2] = MGNFWI(reshape(m0,n),Dobs,wavelet,model,opts,2);
save results_2 results_2


%%
open MGNFWI_camenbert
open GN_up
open pqnl1_2
open minconf_pqn_2
open lbfgsupdate