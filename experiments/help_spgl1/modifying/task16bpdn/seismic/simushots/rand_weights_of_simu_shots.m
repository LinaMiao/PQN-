% random weights of simu shots

clear;close all
%% addpath for PQN working
cd ../../../../../../functions;
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task16bpdn/seismic/simushots/
cd ./simu_functions/
addpath(genpath(pwd))
cd ..

%% set rand generator and sedd
rng(3216,'twister');
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
D = D(1:end,30,1:end);

% Define new data sizes
[nt,nr,ns] = size(D);

% Vectorize D
D = D(:);

% Display
figure
imagesc(reshape(D,nt,ns)); colormap(gray); colorbar;
title('Original data (receiver gather)');
xlabel('Shot number'); ylabel('Time sample')



%% random weights of simu shots
p = .8;
nse = round(p*ns);
% SS = rand(nse,ns);
% Dt = opDirac(nt);
% Dr = opDirac(nr);
% RM = opKron(SS,Dr,Dt);
RM = opKron(opGaussian(ns, nse)',opDirac(nt));

x_test = rand(size(RM,2),1);
y_test = rand(size(RM,1),1);
left = y_test'*(RM*x_test);
right = (RM'*y_test)'*x_test;
error = norm(left-right);
fprintf('In dottest error:%5.5e\n',error);


simD1 = RM*D;
figure;
imagesc(reshape(simD1,nt,nse)); colormap(gray); colorbar;


%% sparsifying transform
% Use this to create a Curvelet SPOT operator:
%C = opCurvelet(nt, ns);
C = opDFT2(nt,ns);

% Transform the data into the Curvelet domain and plot the sorted coefficients 
C_D = C*D;
sort_CD = sort(abs(C_D),'descend');
figure;plot(sort_CD);title('sorted curvelet coefficients')


%% reconstruct
A = RM*C';

% BPDN
options = spgSetParms('optTol', 1e-4, 'iterations', 200);%, 'fid', fid); 
options.fid = fopen('rand_weights_shots_spg.txt','w');
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, simD1, 0, 0, zeros(size(A,2),1), options); % Find BP sol'n.
options.fid = fopen('rand_weights_shots_pqn.txt','w');
sigma_ref = info_spg.rNorm * .99;
[x_pqn1,r_pqn1,g_pqn1,info_pqn1] = pqnl1_2(A, simD1, 0, 0, zeros(size(A,2),1), options,sigma_ref); % Find BP sol'n.
fspg = C'*x_spg;
snrspg = SNR(D,fspg);
fpqn = C'*x_pqn1;
snrpqn = SNR(D,fpqn);

h = figure; 
subplot(1,2,1);imagesc(real(reshape(fspg,nt,ns))); colormap(gray);
title(strcat(['p = .5, SNR=' num2str(snrspg) 'dB']))
subplot(1,2,2);imagesc(real(reshape(fspg-D,nt,ns))); colormap(gray);
title('difference')
saveas(h,'rand_weights_shots_spg.jpg')


h = figure; 
subplot(1,2,1);imagesc(real(reshape(fpqn,nt,ns))); colormap(gray);
title(strcat(['p = .5, SNR=' num2str(snrpqn) 'dB']))
subplot(1,2,2);imagesc(real(reshape(fpqn-D,nt,ns))); colormap(gray);
title('difference')
saveas(h,'rand_weights_shots_pqn.jpg')

% show result
h = figure('Name','Solution paths');
plot(info_spg.xNorm1,info_spg.rNorm2,info_pqn1.xNorm1,info_pqn1.rNorm2);hold on
scatter(info_spg.xNorm1,info_spg.rNorm2);
scatter(info_pqn1.xNorm1,info_pqn1.rNorm2);hold off
legend('SPGL1','PQNl1')
axis tight
saveas(h,'solution_path_rand_weights.jpg');

