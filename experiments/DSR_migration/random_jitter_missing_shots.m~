%Sequential-source data reconstruction (acquistion with randomly jittered missing shots)
clear; close all;
cd ./functions
addpath(genpath(pwd))
cd ..
cd ./simu_functions/
addpath(genpath(pwd))
cd ..

%% original data
% Number of time samples
nt = 1024;
% Number of sources
ns = 178;
% Number of receivers
nr = 178;

% Read data
D = ReadSuFast('GulfOfSuez178.su');
D = reshape(D,nt,nr,ns);

% Select small subset
D = D(1:256,30,1:50);

% Define new data sizes
[nt,nr,ns] = size(D);

% Vectorize D
D = D(:);

% Display
figure
imagesc(reshape(D,nt,ns)); colormap(gray); colorbar;
title('Original data (receiver gather)');
xlabel('Shot number'); ylabel('Time sample')



%% random jittering missing shots
% random jittering missing shots
n = ns;
p = .5;
I_jitter = jitter1d(n,p*n);
S_jitter = zeros(n,1); S_jitter(I_jitter) = 1;
Js = opDiag(S_jitter);
Dt = opDirac(nt);
Dr = opDirac(nr);
RM = opKron(Js,Dr,Dt);


x_test = rand(size(RM,2),1);
y_test = rand(size(RM,1),1);
left = y_test'*(RM*x_test);
right = (RM'*y_test)'*x_test;
error = norm(left-right);
fprintf('In dottest error:%5.5e\n',error);


simD1 = RM*D;
figure;
imagesc(reshape(simD1,nt,ns)); colormap(gray); colorbar;


%% sparsifying transform
% Use this to create a Curvelet SPOT operator:
C = opCurvelet(nt, ns);

% Transform the data into the Curvelet domain and plot the sorted coefficients 
C_D = C*D;
sort_CD = sort(abs(C_D),'descend');
figure;plot(sort_CD);title('sorted curvelet coefficients')


%% reconstruct
options = spgSetParms('optTol', 1e-4, 'iterations', 200);%, 'fid', fid); 
A = RM*C';
xestspg = spgl1_origin(A,simD1,0,1e-3,[],options);
options.iterations = 100;
xestpqn = pqnl1_2(A,simD1,0,1e-3,[],options);
fspg = C'*xestspg;
snrspg = SNR(D,fspg);
fpqn = C'*xestpqn;
snrpqn = SNR(D,fpqn);


    
figure; 
subplot(1,2,1);imagesc(reshape(fspg,nt,ns)); colormap(gray);
title(strcat(['p = .5, SNR=' num2str(snrspg) 'dB']))
subplot(1,2,2);imagesc(reshape(fspg-D,nt,ns)); colormap(gray);
title('difference')


figure; 
subplot(1,2,1);imagesc(reshape(fpqn,nt,ns)); colormap(gray);
title(strcat(['p = .5, SNR=' num2str(snrpqn) 'dB']))
subplot(1,2,2);imagesc(reshape(fpqn-D,nt,ns)); colormap(gray);
title('difference')
