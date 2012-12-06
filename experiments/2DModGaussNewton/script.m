%%
cd ../../../../functions/
addpath(genpath(pwd));
cd ../experiments/2DModGaussNewton
addpath(genpath(pwd));

addpath(genpath('/users/installs/slic/mat_toolbox/spot-slim/'));
addpath(genpath('/users/slic/linamiao/Documents/Documents/Tools/Matlabtools/tristan/'))
%MGNFWI_camenbert

%% inversion
clear
load camenbert
opts.iterations = 50;
opts.fid = fopen('spg','w');
[results_1,~,~,~,infospg] = MGNFWI(reshape(m0,n),Dobs,wavelet,model,opts,1);
save results_1 results_1
save info infospg
clear
load camenbert
load info
sigma_ref = infospg.rNorm;
opts.iterations = 50;
opts.fid = fopen('pqn','w');
[results_2,~,~,~,infopqn] = MGNFWI(reshape(m0,n),Dobs,wavelet,model,opts,2,sigma_ref);
save results_2 results_2


%%
open MGNFWI_camenbert
open GN_up
open lbfgsupdate