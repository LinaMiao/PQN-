%% Sequential-source data reconstruction (acquistion with randomly jittered missing shots)
% The exercise deals with the reconstruction of a fully sampled data-volume
% from data that was subsampled by firing off a few randomly jittered
% shots.
%
% First, we describe the whole process on a receiver gather.

%% addpath
cd ../../../../../../../functions;
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task16bpdn/seismic/simushots/
addpath(pwd)
cd ./simu_functions/
addpath(genpath(pwd))
cd ../haneet


%% Data
%
% Load the data |GulfOfSuez178.su|. We won't read the headers this time and
% set the dimensions manually.
 
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
D = squeeze(D(1:256,30,1:50));
 
% Define new data sizes
[nt,ns] = size(D);    

% Vectorize D
D = D(:);
 
% Display
figure
imagesc(reshape(D,nt,ns)); colormap(gray); colorbar;
title('Original data (receiver gather)'); 
xlabel('Shot number'); ylabel('Time sample')                 


%% Set the parameters for random jittered experiment
% We now subsample the receiver gather by reducing the number of shots
% fired. One such experiment is created by firing the shots at randomly
% jittered locations.

% * Construct the sampling operator |RM| for |p = 0.5| that works on the vectorized version of the data using |opKron| and |opRestriction|. 
% HINT: |opDirac| acts along the time dimension. Use the Kronecker property for the matrix product given in Lab 6.  
% * Test the sampling operator with the dottest.
% * Generate randomly jittered data |jittD| and display the result.
% * Use the adjoint (of the sampling operator) as inverse and display the result. What do you observe?
% * Repeat for different subsampling ratio's.

% Subsampling ratio
p = .8;                        

% Number of randomly jittered shots
nsjitt = round(p * ns);      

% Sampling operator    
idx      = jitter1d(50, 1);   % jitter1d(width, spacing)
idx_perm = randperm(length(idx));
idx      = sort(idx(idx_perm(1:nsjitt)));
RM       = opKron(opRestriction(ns, idx), opDirac(nt)); 

% Randomly jittered data
jittD = RM*D;

% Plot
figure
imagesc(reshape(jittD, nt, nsjitt)); colormap(gray); colorbar;
title('Randomly jittered data'); 
xlabel('Shot number'); ylabel('Time sample')
 
% Adjoint as inverse
figure
imagesc(reshape(RM'*jittD,nt,ns)); colormap(gray); colorbar;
title('Adjoint as inverse')
xlabel('Shot number'); ylabel('Time sample')


%% Sparsifying basis
% To be able to reconstruct the original data from the subsampled data, we
% need a basis in which the data is sparse. We use Curvelets for this. 
% We can think of curvelets as being a localized fourier transform.
%
% First, lets look at the curvelet transform. We can define the curvelet
% transform operator for an image of size |nt| x |ns| as follows:
 
% Use this to create a Curvelet SPOT operator:
 
C = opCurvelet(nt, ns);
 

% %%
% % * Plot a few rows of |C| to get an idea of what a Curvelet looks like.
% % 
%  
% % Define a few unit vectors
% e1 = zeros(size(C,1),1); e1(2000) = 1;
% e2 = zeros(size(C,1),1); e2(20000) = 1;
%  
% % Curvelet in the spatial domain
% figure;
% imagesc(reshape(C'*e1, nt, ns)); colormap(gray)
% title('Curvelet')
%  
% figure;
% imagesc(reshape(C'*e2, nt, ns)); colormap(gray)
% title('Curvelet')


%%
% % Next, we look at how sparse the data is actually is in the Curvelet
% % domain. We can do this by looking at the magnitude of the curvelet coefficients of the
% % data. 
% %
% % * Transform the data into the Curvelet domain and plot the sorted
% % coefficients (use |sort| to sort the coefficients in decreasing order of absolute value).
%  
% % Curvelet coefficients
% x = C*D;
%  
% % Sort the coefficients in descending order, and display
% x_sort = sort(abs(x), 'descend'); 
% figure
% plot(x_sort); axis tight; 
% title('Sorted curvelet coefficients'); xlabel('Index'); ylabel('abs(x)')
 

%% Exercises
% Now, reconstruct the orinigal data from the subsampled data |jittD| for
% different subsampling ratio's:
%
% * Construct the measurement operator |A|.
% HINT: See 'Constructing a suitable matrix' in Lab 7.     
% * Using |spgl1|, estimate the curvelet coefficients |xest|. You can set the 'options' as follows: 
% |options = spgSetParms('optTol', 1e-4, 'iterations', 200, 'fid', fid)|;
% where |fid = fopen('log.txt', 'w');| is the file identifier given to the log file.
% * Recover the data , compute the signal-to-noise ratio.                
% * Display the results (recovered data and the difference).

% Measurement matrix
A = RM*C';
 
% Set spgl1 options
fid     = fopen('log.txt','w');
options = spgSetParms('optTol', 1e-4, 'iterations', 1200, 'fid', fid);
 
% Run spgl1
xest = spgl1(A,jittD,0,0,[],options);
 
% Transform back to physical domain
Dest = C'*xest;
 
% SNR
SNR = -20*log10(norm(D-Dest)/norm(D));
 
% Plot results
figure
imagesc(reshape(Dest,nt,ns)); colormap(gray); colorbar;
title(['Recovered data, SNR = ' num2str(SNR)]); 
xlabel('Shot number'); ylabel('Time sample');
 
figure
imagesc(reshape(Dest - D,nt,ns)); colormap(gray); colorbar;
title('Difference'); 
xlabel('Shot number'); ylabel('Time sample');


%% The bigger picture
% Now lets do things in |3D|.
% PROJECT: Repeat the experiment using the entire data cube.
% Working with |3D| data requires the setup of the sampling operator |RM| and the
% measurement operator |A| to be modified.
% |RM| now incorporates |opDirac| acting along time and receiver dimensions.
% For the measurement operator |A|, the sparse basis |S| is now a
% Kronecker product between the |2D| curvelet (|opCurvelet|) acting across the
% receiver-source plane, and the |1D| wavelet operator (|opSplineWavelet|)
% acting along the time dimension. This Kronecker product is actually a
% parallel Kronecker product (|oppKron2Lo|) between curvelets and wavelets,
% which means parallel computations are performed on multiple processors.
% The operator |S| is setup as:
opW = opSplineWavelet(nt, 1, nt, 3, 5); 

% Determine size of the wavelet operator
nm = opW([],0);      

% Return a SPOT operator using constructor opFunction
W = opFunction(nm{1}, nm{2}, opW);                   

% oppKron2Lo : kronecker tensor product to act on a distributed vector
S = oppKron2Lo(C, W', 1); 


%% Instructions for using the cluster
% When working with the entire data cube, run the experiment in batch mode.
% Do NOT run it interactively. SLIM members can see the
% 'batchsubmit_example.m' script in the
% /users/installs/slic/matlab_etc/contrib folder for details on how to
% submit a batch job. For data of size |256 x 100 x 100 (nt x nr x ns)|,
% use the 'torque1x2' configuration. On average, the batch job runs for 120
% minutes with 200 |spgl1| iterations. Note: do NOT run |spgl1| to
% completion, i.e., fix the number of iterations to get the results in
% reasonable amount of time.


