%Sequential-source data reconstruction (acquistion with randomly jittered missing shots)
clear;close all
%% fix rand seeds
rng(3216,'twister');
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
D = squeeze(D(1:512,30,1:60));

% Define new data sizes
%[nt,nr,ns] = size(D);
[nt,ns] = size(D); 

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
p = .75;
% I_jitter = jitter1d(n,p*n);
% S_jitter = zeros(n,1); S_jitter(I_jitter) = 1;
% Js = opDiag(S_jitter);
% Dt = opDirac(nt);
% Dr = opDirac(nr);
% RM = opKron(Js,Dr,Dt);
nsjitt = round(p * ns);      
idx      = jitter1d(n, n*p);   % jitter1d(width, spacing)
idx_perm = randperm(length(idx));
idx      = sort(idx(idx_perm(1:nsjitt)));
RM       = opKron(opRestriction(ns, idx), opDirac(nt)); 


x_test = rand(size(RM,2),1);
y_test = rand(size(RM,1),1);
left = y_test'*(RM*x_test);
right = (RM'*y_test)'*x_test;
error = norm(left-right);
fprintf('In dottest error:%5.5e\n',error);


simD1 = RM*D;
h = figure;
imagesc(reshape(simD1,nt,round(p*ns))); colormap(gray); colorbar;
saveas(h,'jitter_missing_receivergather')


%% sparsifying transform
% Use this to create a Curvelet SPOT operator:
C = opCurvelet(nt, ns);

% Transform the data into the Curvelet domain and plot the sorted coefficients 
C_D = C*D;
sort_CD = sort(abs(C_D),'descend');
figure;plot(sort_CD);title('sorted curvelet coefficients')


%% reconstruct
options = spgSetParms('optTol', 1e-4, 'iterations', 500);%, 'fid', fid); 
A = RM*C';

% BPDN
options.fid = fopen('jitter_missing_spg.txt','w');
t = cputime; 
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, simD1, 0, 0, zeros(size(A,2),1), options); % Find BP sol'n.
tspg = cputime - t;
fspg = C'*x_spg;
snrspg = SNR(D,fspg);
save infospg x_spg tspg info_spg fspg snrspg
options.fid = fopen('jitter_missing_pqn.txt','w');
sigma_ref = info_spg.rNorm;

corrections_list = [5 10 15 20 25 30];
for i = 1:length(corrections_list)
    corrections = corrections_list(i);
    t = cputime;
    [x_pqn1(:,i),r_pqn1,g_pqn1,info_pqn1(:,i)] = pqnl1_2(A, simD1, 0, 0, zeros(size(A,2),1), options,sigma_ref,corrections); % Find BP sol'n.
    tpqn(:,i) = cputime - t;
    
    
    fpqn(:,i) = C'*x_pqn1(:,i);
    snrpqn(:,i) = SNR(D,fpqn(:,i));
    
    save info_pqn x_pqn1 tpqn info_pqn1 fpqn snrpqn
end

h = figure; 
subplot(1,2,1);imagesc(reshape(fspg,nt,ns)); colormap(gray);caxis([-400 400])
title(strcat(['p = .75, SNR=' num2str(snrspg) 'dB']))
subplot(1,2,2);imagesc(reshape(fspg-D,nt,ns)); colormap(gray);caxis([-400 400])
title('difference')
saveas(h,'jitter_missing_spg')

h = figure; 
subplot(1,2,1);imagesc(reshape(fpqn(:,i),nt,ns)); colormap(gray);caxis([-400 400])
title(strcat(['p = .75, SNR=' num2str(snrpqn(:,i)) 'dB']))
subplot(1,2,2);imagesc(reshape(fpqn(:,i)-D,nt,ns)); colormap(gray);caxis([-400 400])
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

