%% addpath of tools and functions
clear; close all;
cd ./simu_functions/
addpath(genpath(pwd))
cd ../../../../../../functions/;
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task11lasso/seismic/


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

% calculate tau
options = spgSetParms('optTol', 1e-4, 'iterations', 1000);%, 'fid', fid); 
xestspg = spgl1(A,simD1,0,1e-3,[],options);
tau = norm(xestspg,1);
%tau = 2.2072179e+05;



% spg and pqn
options = spgSetParms('optTol', 1e-4, 'iterations', 200);%, 'fid', fid); 
options.fid = fopen('jitter_missing_spg.txt','w');
xinit = zeros(size(A,2),1);

[xestspg,~,~,infospg] = spgl1(A,simD1,tau,[],xinit,options);
options.fid = fopen('jitter_missing_pqn.txt','w');
options.optTol = infospg.rNorm;
xestpqn = pqnl1_2(A,simD1,tau,[],xinit,options);
fspg = C'*xestspg;
snrspg = SNR(D,fspg);
fpqn = C'*xestpqn;
snrpqn = SNR(D,fpqn);



    
h = figure; 
subplot(1,2,1);imagesc(reshape(fspg,nt,ns)); colormap(gray);
title(strcat(['p = .5, SNR=' num2str(snrspg(1)) 'dB']))
subplot(1,2,2);imagesc(reshape(fspg-D,nt,ns)); colormap(gray);
title('difference')
saveas(h,'jitter_missing_spg.jpg')

h = figure; 
subplot(1,2,1);imagesc(reshape(fpqn,nt,ns)); colormap(gray);
title(strcat(['p = .5, SNR=' num2str(snrpqn(1)) 'dB']))
subplot(1,2,2);imagesc(reshape(fpqn-D,nt,ns)); colormap(gray);
title('difference')
saveas(h,'jitter_missing_pqn.jpg');

%% if given known strict sparse vector
[m n] = size(A); k = .2*round(n/log(m));
p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
figure;plot(x0)
b0  = A*x0;

tau = norm(x0,1);

options = spgSetParms('optTol', 1e-4, 'iterations', 200);%, 'fid', fid);
options.fid = fopen('jitter_missing_sparse_spg,txt','w');
xinit = zeros(size(A,2),1);

[xestspg,~,~,infospg] = spgl1(A,b0,tau,[],xinit,options);
tau = infospg.rNorm;
options.fid = fopen('jitter_missing_sparse_pqn,txt','w');
xestpqn = pqnl1_2(A,b0,tau,[],xinit,options);
snrspg = SNR(x0,xestspg);
snrpqn = SNR(x0,xestpqn);

h = figure('Name','strcit sparse vector SPG'); 
subplot(2,1,1);plot(xestspg); 
title(strcat(['p = .5, SNR=' num2str(snrspg) 'dB']))
subplot(2,1,2);plot(xestspg - x0);
title('difference')
saveas(h,'jitter_missing_sparse_spg,jpg');

h = figure('Name','strcit sparse vector PQN'); 
subplot(2,1,1);plot(xestpqn); 
title(strcat(['p = .5, SNR=' num2str(snrpqn) 'dB']))
subplot(2,1,2);plot(xestpqn - x0);
title('difference')
saveas(h,'jitter_missing_sparse_pqn.jpg');

%% if given known compressible vector
nn = linspace(0,1,n);
x0_compress = exp(-nn.^.1);
x0_compress = x0_compress - min(x0_compress);
figure;plot(x0_compress)
x0_compress = x0_compress(:);
b_compress  = A*x0_compress + 0.005 * randn(m,1);


tau = norm(x0_compress,1);

options = spgSetParms('optTol', 1e-4, 'iterations', 200);%, 'fid', fid); 
options.fid = fopen('jitter_missing_compress_spg','w');
xinit = zeros(size(A,2),1);

[xestspg,~,~,infospg] = spgl1(A,b_compress,tau,[],xinit,options);
tau = infospg.rNorm;
options.fid = fopen('jitter_missing_compress_pqn','w');

xestpqn = pqnl1_2(A,b_compress,tau,[],xinit,options);
snrspg = SNR(x0_compress,xestspg);
snrpqn = SNR(x0_compress,xestpqn);

h = figure('Name','compressible vector SPG'); 
subplot(2,1,1);plot(xestspg); 
title(strcat(['p = .5, SNR=' num2str(snrspg) 'dB']))
subplot(2,1,2);plot(xestspg - x0_compress);
title('difference')
saveas(h,'jitter_missing_compress_spg.jpg');

h = figure('Name','compressible vector PQN'); 
subplot(2,1,1);plot(xestpqn); 
title(strcat(['p = .5, SNR=' num2str(snrpqn) 'dB']))
subplot(2,1,2);plot(xestpqn - x0_compress);
title('difference')
saveas(h,'jitter_missing_compress_pqn,jpg');


