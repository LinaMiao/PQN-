%%
cd ../../../../functions/
addpath(genpath(pwd));
cd ../experiments/2DModGaussNewton/
addpath(genpath(pwd));
cd ./examples/example_BG/
addpath(genpath('/users/installs/slic/mat_toolbox/spot-slim/'));
addpath(genpath('/users/slic/linamiao/Documents/Documents/Tools/Matlabtools/tristan/'))



% This example produce the result from Modified Gauss-Newten Full-waveform
% inversion.
% -----------------------------------------------
%
% In this example, we process FWI from low frequency band to high frequency
% band. 10 GN subproblems are solved for each frequency band which contains 
% 10 frequencies. For each GN subproblem, we use 7 randomly selected simultaneous
% shots and roughly 20 L1 solver iterations
% 
% We sugguest you run this script with parallel matlab, 10 workers will be the best
%
% Results will be saved as a mat file under results directory.
% ------------------------------------------------
%
% setup a random seed

%% setup parameters for the experiment
clear;  close all;
model.minf        = 2.9; % minimal frequency
model.maxf        = 22;  % maximal frequency
model.sourcet     = 1; % source type: 1 sim source; 2 sequential source.
model.name        = ['/example_BG/FWIresult_BG']; % name of result
model.vmin  = 1480;model.vmax  = 4800;% min vel and max vel for the projection
model.nf          = 10; % number of frequencies in each frequency band
model.snf         = 10; % number of frequencies for each GN updates 
model.ol          = 5;  % number of overlap frequencies for two adjacent bands
model.gl          = 10; % grid length
model.nshot       = 351;% number of shot positions
model.nsim        = 7; % Number of sim-shots for each GN subproblem
model.sdep        = 3; % source depth, unit: meters
model.sp          = round(linspace(1,701,351)); % source position
model.rdep        = 2; % receiver depth, unit: meters
model.nrec        = 701; % number of receivers
model.rp          = round(linspace(1,701,701)); % receiver positions
model.water       = 16; % estimated water depth, unit: number of grid
save model model
%% load observation data and wavelet information
load('obdata.mat')
Dobs        = data; clear data % observation data 
wavelet     = wavelet;         % wavelet info



% load initial model and wavelet info

%% inversion

opts.iterations = 50;
opts.fid = fopen('spg','w');
load bg2dmodel.mat
t = cputime;
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
[results_1,~,~,~,infospg] = MGNFWI(vel1,Dobs,wavelet,model,opts,1);
t_spg = cputime - t;
save results_1 results_1
save infospg infospg

% pqn
clear
load model
load('obdata.mat')
Dobs        = data; clear data % observation data 
wavelet     = wavelet;         % wavelet info
load bg2dmodel.mat
% sigma_ref
load infospg
sigma_ref = zeros(50,1);
for i = 1:length(infospg)
    sigma_ref(i) = infospg{i}.rNorm;
end
opts.iterations = 50;
opts.fid = fopen('pqn','w');
t = cputime;
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
[results_2,~,~,~,infopqn] = MGNFWI(vel1,Dobs,wavelet,model,opts,2,sigma_ref);
t_pqn = cputime - t;
save results_2 results_2
save infopqn infopqn