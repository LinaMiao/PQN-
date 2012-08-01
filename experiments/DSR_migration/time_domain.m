%% 
% time domain.

%addpath of some functions and operators
% addpath /Volumes/Users/linamiao/Documents/Tools/Matlabtools/tuning/
addpath(genpath('./simu_functions'))
addpath ../data
addpath(genpath('../functions'))

clear
% close all
%% read data
[data,SuTraceHeaders,SuHeader]=ReadSu('data_ex3.su');

% source and receiver coordinates per trace
for k=1:length(SuTraceHeaders)
    xs(k)= SuTraceHeaders(k).SourceX;
    xr(k)= SuTraceHeaders(k).GroupX;
end

% note that the sample interval is given in miliseconds!
t = [0:SuHeader.ns-1]'*SuHeader.dt*1e-6;

% source and receiver coordinate vectors
xs = unique(xs); 
xr = unique(xr); 

% reshape data into cube
data = reshape(data,length(t),length(xr),length(xs));

% window and subsample
data(1:100,:,:) = 0;
data = data(1:end,1:10:end,:);
data = data(:,1:end-1,1:end-1);%SR2MH cannot handle odd no.
xs = xs(1:end-1); 
xr = xs;
t = t(1:end);

%data(1:100,:,:) = 0;
%data = data(:,1:10:end,:);
%xr = xs;

nt = length(t);
ns = length(xs);nsrc = length(xs);
nr = length(xr);nrec = length(xr);


%% simulate simutaneous experiment
% % measurement matrix
% p = .5;
% RM = opSimSourcePeriodTimeDither([nt nr ns], [floor(p*nt*ns),nr], ns);
% 
% % observed data from simutaneous experiment
% simD = RM*vec(data); 
% simd = reshape(simD,64,11,11);
% t = t(1:nt/2);

%% simutaneous experiment with randomized source superposition
% source superposition
ncsrc = 5;
M   = opGaussian(nsrc,ncsrc);
RM  = opKron(opDirac(nt),opDirac(nrec),M');

simD = RM*vec(data);
simD = reshape(simD,nt,nrec,ncsrc);
%% sparsity transform in wavefiled domain
% % wavelet operator along time axis
% W = opSplineWavelet(nt, 1, nt, 3, 5);
% 
% % curvelet operator along source-receiver oordinates
% C = opCurvelet(nr,ns,6,16,1,'ME',0);
% 
% % oppKron2Lo : kronecker tensor product to act on a distributed vector
% S = opKron(C, W', 1); 

%% DSR extended migration for regular sources
z = 10:10:1000; nz = length(z);
v = 2000*ones(length(z),1);

dim = [nt,nrec,nsrc];
K = opDSR(dim,t,xr,xs,z,v);
M = K*vec(data); 


% % check1
% a = test11(data,t,xr,xs,z,v(1:length(z)));
% image = zeros(nz,nrec);
% for i = 1:100
%     image(i,:) = diag(squeeze(a(i,:,:)));
% end
% figure;imagesc(real(image));
% 
% 
% 
% % check 2
% [output] = imgather(a,z,xr,xs,t); % output all midpoint imagather
% nm = size(output,1);
% 
% figure;
% for i = 1:nm
%     imageHandle = subplot(1,nm,i); imagesc(real(squeeze(output(i,:,:))));colormap(gray);axis off;
%     pHandle = get(imageHandle, 'pos');
%     pHandle(3) = pHandle(3) + 0.005;
%     set(imageHandle, 'pos', pHandle);
% end


% check
ny = nr; % number of midpoints
np = numel(M)/(ny*nz); % number of wavenumbers
M = reshape(M,ny,nz,np); 
figure;
for i = 1:ny
    imageHandle = subplot(1,ny,i); imagesc(real(squeeze(M(i,:,:))));colormap(gray);axis off;
    pHandle = get(imageHandle, 'pos');
    pHandle(3) = pHandle(3) + 0.005;
    set(imageHandle, 'pos', pHandle);
end
% colorbar;





%% DSR extended migration for sim experiments
% % DSR result if sim
% % periodical locate simu sources, the corresponding source grid
% % xcs = linspace(xs(1),xs(end),ns);
% 
% % apply the adjoint of simu experiment matrix
% Data1 = RM'*vec(simD);
% 
% % perform DSR 
% M1 = K*Data1; 
% 
% % check
% %M1 = reshape(M1,ny,nz,np); 
% M1 = reshape(M1,ny*nz,np);
% figure;imagesc(M1);colormap(gray);
% 
% 
%% sparsity transform in image domain
% in dirac basis
S = opKron(opDirac(np),opDirac(nz*ny));
% norm(M(:,1),1)
% norm(M1(:,1),1)
% % compute the ||X||_1,2 of M and M1
% Norm12M = mixNorm(M,1,2)
% Norm12M1 = mixNorm(M1,1,2)
% 
% 
% % in fourier basis
% norm(fft(M(:,1)),1)
% norm(fft(M1(:,1)),1)
% % compute the ||X||_1,2 of fftM and fftM1
% Norm12M = mixNorm(fft(M,[],1),1,2)
% Norm12M1 = mixNorm(fft(M1,[],1),1,2)
% % 
%% inversion via spgl1
A = RM*K'*S';
B = vec(simD);
fid = fopen('log_log.txt', 'w'); 
opts = spgSetParms('optTol',1e-4,'iterations',20);
sigma = 1e-3;
A_handle = @(x,mode)A_handle_multiply(A,x,mode);
[x_hat,~,~,info] = spgl1(A,B,0,sigma,[],opts);

%% inversion via joint sparsity
% %[X_hat,R,G,INFO] = spg_mmv(A_handle,B,sigma,opts);
% 
% 
%% evaluate result
% m1 = S'*x_hat;
% exImage = K'*vec(m1);
% op = opTEST11(dim,t,xr,xs,z,v);
% exImage = op*exImage;
% exImage = reshape(exImage,nz,nrec,nsrc);
% 
% image = ieximage(exImage);
% figure;imagesc(image);
% 
% 
% 
% 
% 
% 
% M1 = S'*X_hat;
% SNR1 = snr(vec(M),vec(M1));
% 
% D1 = K'*M1;
% SNR2 = snr(vec(D1),vec(data));
% 
