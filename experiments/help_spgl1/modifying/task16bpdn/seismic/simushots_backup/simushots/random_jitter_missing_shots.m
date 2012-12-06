%Sequential-source data reconstruction (acquistion with randomly jittered missing shots)
clear;close all
%% addpath for PQN working
cd ../../../../../../functions;
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task16bpdn/seismic/simushots/
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
D = D(1:256,30,:);

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
p = .6;
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
h = figure;
imagesc(reshape(simD1,nt,ns)); colormap(gray); colorbar;
saveas(h,'jitter_missing_receivergather')


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

% BPDN
options.fid = fopen('jitter_missing_spg.txt','w');
t = cputime; 
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, simD1, 0, 0, zeros(size(A,2),1), options); % Find BP sol'n.
tspg = cputime - t;
save infospg x_spg tspg info_spg 
options.fid = fopen('jitter_missing_pqn.txt','w');
sigma_ref = info_spg.rNorm;
t = cputime;
[x_pqn1,r_pqn1,g_pqn1,info_pqn1] = pqnl1_2(A, simD1, 0, 0, zeros(size(A,2),1), options,sigma_ref); % Find BP sol'n.
tpqn = cputime - t;
save info_pqn x_pqn1 tpqn info_pqn1
fspg = C'*x_spg;
snrspg = SNR(D,fspg);
fpqn = C'*x_pqn1;
snrpqn = SNR(D,fpqn);

h = figure; 
subplot(1,2,1);imagesc(reshape(fspg,nt,ns)); colormap(gray);
title(strcat(['p = .6, SNR=' num2str(snrspg(1)) 'dB']))
subplot(1,2,2);imagesc(reshape(fspg-D,nt,ns)); colormap(gray);
title('difference')
saveas(h,'jitter_missing_spg')

h = figure; 
subplot(1,2,1);imagesc(reshape(fpqn,nt,ns)); colormap(gray);
title(strcat(['p = .6, SNR=' num2str(snrpqn(1)) 'dB']))
subplot(1,2,2);imagesc(reshape(fpqn-D,nt,ns)); colormap(gray);
title('difference')
saveas(h,'jitter_missing_pqn');

% show result
h = figure('Name','Solution paths');
plot(info_spg.xNorm1,info_spg.rNorm2,info_pqn1.xNorm1,info_pqn1.rNorm2);hold on
scatter(info_spg.xNorm1,info_spg.rNorm2);
scatter(info_pqn1.xNorm1,info_pqn1.rNorm2);hold off
legend('SPGL1','PQNl1')
axis tight
saveas(h,'solution_path_jitter_missing');

