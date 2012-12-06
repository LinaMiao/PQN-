    
%% Installation
% download and install 
clear; close all;
cd ./simu_functions/
addpath(genpath(pwd))
cd ../..
addpath(genpath(pwd))
cd ../../../../pqnl1;
addpath(genpath(pwd))

cd ../experiments/help_spgl1/modifying/task11lasso/seismic
rmpath('/Volumes/Users/linamiao/Dropbox/PQN/pqnl1/minConF/')
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
p_list = [.5];

p = p_list;
D_RM1 = opSimSourceRandTimeDither([nt,nr,ns],[p*nt*ns,1],ns);
simD1 = D_RM1*D;
A = D_RM1*C';

options = spgSetParms('optTol', 1e-4, 'iterations', 1000);%, 'fid', fid); 
xestspg = spgl1(A,simD1,0,1e-3,[],options);
tau = norm(xestspg,1);
%tau = 2.2072179e+05;

options = spgSetParms('optTol', 1e-4, 'iterations', 200);%, 'fid', fid); 
xinit = zeros(size(A,2),1);

which spgl1
%keyboard;
xestspg = spgl1(A,simD1,tau,[],xinit,options);
%options.iterations = 100;
xestpqn = pqnl1_2(A,simD1,tau,[],xinit,options);
fspg = C'*xestspg;
snrspg = SNR(D,fspg);
fpqn = C'*xestpqn;
snrpqn = SNR(D,fpqn);



    
figure; 
subplot(1,2,1);imagesc(reshape(fspg,nt,ns)); colormap(gray);
title(strcat(['p = .5, SNR=' num2str(snrspg(1)) 'dB']))
subplot(1,2,2);imagesc(reshape(fspg-D,nt,ns)); colormap(gray);
title('difference')


figure; 
subplot(1,2,1);imagesc(reshape(fpqn,nt,ns)); colormap(gray);
title(strcat(['p = .5, SNR=' num2str(snrpqn(1)) 'dB']))
subplot(1,2,2);imagesc(reshape(fpqn-D,nt,ns)); colormap(gray);
title('difference')


%% if given known strict sparse vector
[m n] = size(A); k = .2*round(n/log(m));
p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
figure;plot(x0)
b0  = A*x0;

tau = norm(x0,1);

options = spgSetParms('optTol', 1e-4, 'iterations', 200);%, 'fid', fid); 
xinit = zeros(size(A,2),1);

xestspg = spgl1(A,b0,tau,[],xinit,options);
xestpqn = pqnl1_2(A,b0,tau,[],xinit,options);
snrspg = SNR(x0,xestspg);
snrpqn = SNR(x0,xestpqn);

figure('Name','strcit sparse vector SPG'); 
subplot(2,1,1);plot(xestspg); 
title(strcat(['p = .5, SNR=' num2str(snrspg) 'dB']))
subplot(2,1,2);plot(xestspg - x0);
title('difference')

figure('Name','strcit sparse vector PQN'); 
subplot(2,1,1);plot(xestpqn); 
title(strcat(['p = .5, SNR=' num2str(snrpqn) 'dB']))
subplot(2,1,2);plot(xestpqn - x0);
title('difference')

%% if given known compressible vector
nn = linspace(0,1,n);
x0_compress = exp(-nn.^.1);
x0_compress = x0_compress - min(x0_compress);
figure;plot(x0_compress)
x0_compress = x0_compress(:);
b_compress  = A*x0_compress + 0.005 * randn(m,1);


tau = norm(x0_compress,1);

options = spgSetParms('optTol', 1e-4, 'iterations', 200);%, 'fid', fid); 
xinit = zeros(size(A,2),1);

xestspg = spgl1(A,b_compress,tau,[],xinit,options);
xestpqn = pqnl1_2(A,b_compress,tau,[],xinit,options);
snrspg = SNR(x0_compress,xestspg);
snrpqn = SNR(x0_compress,xestpqn);

figure('Name','compressible vector SPG'); 
subplot(2,1,1);plot(xestspg); 
title(strcat(['p = .5, SNR=' num2str(snrspg) 'dB']))
subplot(2,1,2);plot(xestspg - x0_compress);
title('difference')

figure('Name','compressible vector PQN'); 
subplot(2,1,1);plot(xestpqn); 
title(strcat(['p = .5, SNR=' num2str(snrpqn) 'dB']))
subplot(2,1,2);plot(xestpqn - x0_compress);
title('difference')


