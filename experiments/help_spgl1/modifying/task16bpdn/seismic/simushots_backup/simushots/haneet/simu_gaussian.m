%% Sequential-source data reconstruction from simultaneous data simulated with a Gaussian matrix
% The exercise deals with the reconstruction of a fully sampled data-volume
% from data that was subsampled by firing off sources simultaneously with
% random weights.
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
D = squeeze(D(1:256,30,:));
 
% Define new data sizes
[nt,ns] = size(D);    
 
% Vectorize D
D = D(:);
 
% Display
figure
imagesc(reshape(D,nt,ns)); colormap(gray); colorbar;
title('Original data (receiver gather)'); 
xlabel('Shot number'); ylabel('Time sample')            
 
 
%% Set the parameters for simultaneous-source experiment
% We now subsample the receiver gather by creating simultaneous source
% experiments. One such experiment is created by multiplying the data
% matrix by a random vector of length |ns|. We can create |nse| of these simultaneous-source
% experiments by multiplying by a random matrix of size |ns| x |nse|. The ration between |nse| and |ns| is
% called the subsampling ratio |p|.
%
% * Construct the sampling operator |RM| for |p = 0.5| that works on the vectorized version of the data using |opKron| and |opGaussian|. 
% HINT: |opDirac| acts along the time dimension. Use the Kronecker property for the matrix product given in Lab 6.  
% * Test the sampling operator with the dottest.
% * Generate simultaneous data |simD| and display the result.
% * Use the adjoint (of the sampling operator) as inverse and display the result. What do you observe?
% * Repeat for different subsampling ratio's.
 
% Subsampling ratio
p = 0.8;                        
 
% Number of simultaneous-source experiments
nse = round(p * ns);      
 
% Sampling operator
RM = opKron(opGaussian(ns, nse)',opDirac(nt));
 
% Simultaneous data
simD = RM*D;
 
% Plot
figure
imagesc(reshape(simD,nt,nse)); colormap(gray); colorbar;
title('Simultaneous data'); 
xlabel('Shot number'); ylabel('Time sample')
 
% Adjoint as inverse
figure
imagesc(reshape(RM'*simD,nt,ns)); colormap(gray); colorbar;
title('Adjoint as inverse'); 
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

 
%%
% * Plot a few rows of |C| to get an idea of what a Curvelet looks like.
% 
 
% Define a few unit vectors
e1 = zeros(size(C,1),1); e1(2000) = 1;
e2 = zeros(size(C,1),1); e2(20000) = 1;
 
% Curvelet in the spatial domain
figure;
imagesc(reshape(C'*e1, nt, ns)); colormap(gray)
title('Curvelet')
 
figure;
imagesc(reshape(C'*e2, nt, ns)); colormap(gray)
title('Curvelet')
 
 
%%
% Next, we look at how sparse the data is actually is in the Curvelet
% domain. We can do this by looking at the magnitude of the curvelet coefficients of the
% data. 
%
% * Transform the data into the Curvelet domain and plot the sorted
% coefficients (use |sort| to sort the coefficients in decreasing order of absolute value).
 
% Curvelet coefficients
x = C*D;
 
% Sort the coefficients in descending order, and display
x_sort = sort(abs(x), 'descend'); 
figure
plot(x_sort); axis tight; 
title('Sorted curvelet coefficients'); xlabel('Index'); ylabel('abs(x)')
 
 
%% Exercises
% Now, reconstruct the orinigal data from the subsampled data |simD| for
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
fid     = fopen('logspg.txt','w');
options = spgSetParms('optTol', 1e-4, 'iterations', 1000, 'fid', fid);
 
% Run spgl1
t = cputime;
[xest,r_spg,g_spg,info_spg] = spgl1(A,simD,0,1e-3,[],options);
tspg = cputime - t;
% Transform back to physical domain
Dest = C'*xest;
 
% SNR
SNRspg = -20*log10(norm(D-Dest)/norm(D));
  
% Plot results
figure
imagesc(reshape(Dest,nt,ns)); colormap(gray); colorbar;
title(['Recovered data, SNR = ' num2str(SNRspg)]); 
xlabel('Shot number'); ylabel('Time sample');
 
figure
imagesc(reshape(Dest - D,nt,ns)); colormap(gray); colorbar;
title('Difference'); 
xlabel('Shot number'); ylabel('Time sample');

 
%% PQN
A = RM*C';
 
% Set spgl1 options
fid     = fopen('logpqn.txt','w');
options = spgSetParms('optTol', 1e-4, 'iterations', 1000, 'fid', fid);
 
% Run spgl1
sigma_ref = info_spg.rNorm;
t = cputime;
[xest,r_pqn1,g_pqn1,info_pqn1] = pqnl1_2(A,simD,0,1e-3,zeros(size(A,2),1),options,sigma_ref);
tpqn = cputime -t ;
% Transform back to physical domain
Dest = C'*xest;
 
% SNR
SNRpqn = -20*log10(norm(D-Dest)/norm(D));
  
% Plot results
figure
imagesc(reshape(Dest,nt,ns)); colormap(gray); colorbar;
title(['Recovered data, SNR = ' num2str(SNRpqn)]); 
xlabel('Shot number'); ylabel('Time sample');
 
figure
imagesc(reshape(Dest - D,nt,ns)); colormap(gray); colorbar;
title('Difference'); 
xlabel('Shot number'); ylabel('Time sample');

save temp

