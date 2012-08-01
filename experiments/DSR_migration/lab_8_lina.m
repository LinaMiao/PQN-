% Sequential-source data reconstruction from randomized 'marine' acquisiton

% this is for the lab(project2) from eosc 454

% lina miao
% 74721119



% The exercise deals with the reconstruction of a fully sampled data-volume
% from data that was subsampled by firing randomly dithered sources in marine.

%% Installation
% download and install 
clear; close all;
cd ./functions
addpath(genpath(pwd))
cd ..
cd ./simu_functions/
addpath(genpath(pwd))
cd ..
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
RM2 = opSimSourcePeriodTimeDither([10 1 10], [5*10 1], 10);

% plot very long time series
figure;
plot(I(:),1:length(I(:)),'*');xlim([0.5 1.5]);ylabel('time samples');

% plot compressed series
figure;
plot(RM1*I(:),1:length(RM1*I(:)),'o');xlim([0.5 1.5]);ylabel('time samples');
figure;
plot(RM2*I(:),1:length(RM2*I(:)),'o');xlim([0.5 1.5]);ylabel('time samples');


%%
% Construct the sampling operator RM for p = 0.5 that works on the vectorized 
% version of the data using opSimSourceRandTimeDither.
p = .5;
D_RM1 = opSimSourceRandTimeDither([nt,nr,ns],[p*nt*ns,1],ns);
D_RM2 = opSimSourcePeriodTimeDither([nt,nr,ns],[p*nt*ns,1],ns);

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
simD2 = D_RM2*D;
figure;
imagesc(reshape(simD1,p*nt,ns)); colormap(gray); colorbar;
figure;
imagesc(reshape(simD2,p*nt,ns)); colormap(gray); colorbar;


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
options = spgSetParms('optTol', 1e-4, 'iterations', 200);%, 'fid', fid); 
p_list = [.5];
for i = 1:1
    p = p_list(i);
    D_RM1 = opSimSourceRandTimeDither([nt,nr,ns],[p*nt*ns,1],ns);
    simD1 = D_RM1*D;
    A = D_RM1*C';
    xestspg(:,i) = spgl1_origin(A,simD1,0,1e-3,[],options);
    options.iterations = 100;
    xestpqn(:,i) = pqnl1_2(A,simD1,0,1e-3,[],options);
    fspg(:,i) = C'*xestspg(:,i);
    snrspg(i) = SNR(D,fspg(:,i)); 
    fpqn(:,i) = C'*xestpqn(:,i);
    snrpqn(i) = SNR(D,fpqn(:,i));
end


    
figure; 
subplot(1,2,1);imagesc(reshape(fspg(:,1),nt,ns)); colormap(gray);
title(strcat(['p = .5, SNR=' num2str(snrspg(1)) 'dB']))
subplot(1,2,2);imagesc(reshape(fspg(:,1)-D,nt,ns)); colormap(gray);
title('difference')


figure; 
subplot(1,2,1);imagesc(reshape(fpqn(:,1),nt,ns)); colormap(gray);
title(strcat(['p = .5, SNR=' num2str(snrpqn(1)) 'dB']))
subplot(1,2,2);imagesc(reshape(fpqn(:,1)-D,nt,ns)); colormap(gray);
title('difference')

