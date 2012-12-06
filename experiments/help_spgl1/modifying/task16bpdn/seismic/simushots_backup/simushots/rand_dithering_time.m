    
clear;close all
%% addpath for PQN working
cd ../../../../../../functions;
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task16bpdn/seismic/simushots/
%% Data
% Number of time samples
nt = 1024;
% Number of sources
ns = 178;
% Number of receivers
nr = 178;

% Time sampling interval
dt = 0.004;

% Read data
D = ReadSuFast('GulfOfSuez178.su');
D = reshape(D,nt,nr,ns);

% Select small subset
D = D(1:256,30,1:100);

% Define new data sizes
[nt,nr,ns] = size(D);

% Vectorize D
D = D(:);

% Display
figure
imagesc(reshape(D,nt,ns)); colormap(gray); colorbar;
title('Original data (receiver gather)');
xlabel('Shot number'); ylabel('Time sample')
%% Set the patameters for randomized experiment
I  = eye(10);
RM1 = opSimSourceRandTimeDither([10 1 10], [5*10 1], 10);

% plot very long time series
figure;
plot(I(:),1:length(I(:)),'*');xlim([0.5 1.5]);ylabel('time samples');

% plot compressed series
figure;
plot(RM1*I(:),1:length(RM1*I(:)),'o');xlim([0.5 1.5]);ylabel('time samples');


%%
% Construct the sampling operator RM for p = 0.5 that works on the vectorized 
% version of the data using opSimSourceRandTimeDither.
p = .5;
D_RM1 = opSimSourceRandTimeDither([nt,nr,ns],[p*nt*ns,1],ns);

% Test the sampling operator with the dottest.
x_test = rand(size(D_RM1,2),1);
y_test = rand(size(D_RM1,1),1);
left = y_test'*(D_RM1*x_test);
right = (D_RM1'*y_test)'*x_test;
error = norm(left-right);
fprintf('In dottest error:%5.5e\n',error);

%%
% Generate simultaneous data simD and display the result.
simD1 = D_RM1*D;
figure;
imagesc(reshape(simD1,p*nt,ns)); colormap(gray); colorbar;


%% sparsifying transform
% Use this to create a Curvelet SPOT operator:
C = opCurvelet(nt, ns);

% Transform the data into the Curvelet domain and plot the sorted coefficients 
C_D = C*D;
sort_CD = sort(abs(C_D),'descend');
figure;plot(sort_CD);title('sorted curvelet coefficients')
%% exercises
% Construct the measurement operator A. HINT: See 'Constructing a suitable 
% matrix' in Lab 7.
% Using spgl1, estimate the curvelet coefficients xest.

fid = fopen('log.txt', 'w'); 
p = [.5];
D_RM1 = opSimSourceRandTimeDither([nt,nr,ns],[p*nt*ns,1],ns);
simD1 = D_RM1*D;
A = D_RM1*C';

% BPDN
options = spgSetParms('optTol', 1e-4, 'iterations', 200);%, 'fid', fid); 
options.fid = fopen('rand_dthering_time_spg.txt','w');
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, simD1, 0, 0, zeros(size(A,2),1), options); % Find BP sol'n.
options.fid = fopen('rand_dithering_time_pqn.txt','w');
sigma_ref = info_spg.rNorm;
[x_pqn1,r_pqn1,g_pqn1,info_pqn1] = pqnl1_2(A, simD1, 0, 0, zeros(size(A,2),1), options,sigma_ref); % Find BP sol'n.
fspg = C'*x_spg;
snrspg = SNR(D,fspg);
fpqn = C'*x_pqn1;
snrpqn = SNR(D,fpqn);

h = figure; 
subplot(1,2,1);imagesc(reshape(fspg,nt,ns)); colormap(gray);
title(strcat(['p = .5, SNR=' num2str(snrspg(1)) 'dB']))
subplot(1,2,2);imagesc(reshape(fspg-D,nt,ns)); colormap(gray);
title('difference')
saveas(h,'rand_dithering_time_spg.jpg')

h = figure; 
subplot(1,2,1);imagesc(reshape(fpqn,nt,ns)); colormap(gray);
title(strcat(['p = .5, SNR=' num2str(snrpqn(1)) 'dB']))
subplot(1,2,2);imagesc(reshape(fpqn-D,nt,ns)); colormap(gray);
title('difference')
saveas(h,'rand_dithering_time_pqn.jpg');

% show result
h = figure('Name','Solution paths');
plot(info_spg.xNorm1,info_spg.rNorm2,info_pqn1.xNorm1,info_pqn1.rNorm2);hold on
scatter(info_spg.xNorm1,info_spg.rNorm2);
scatter(info_pqn1.xNorm1,info_pqn1.rNorm2);hold off
legend('SPGL1','PQNl1')
axis tight
saveas(h,'solution_path_dithering_time.jpg');


