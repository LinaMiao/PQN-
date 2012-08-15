%% addpath for PQN working
cd ../../../../../functions;
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task12illconditioned/convolution/
%stream = RandStream.getGlobalStream;
%reset(stream);


%% problem setting
% time axis
t = [0:.001:2]';
N = length(t);

% true signal g has approx k spikes with random amplitudes
k = 20;
g = zeros(N,1);
g(randi(N,k,1)) = randn(k,1);

% filter
w = (1-2*1e3*(t-.2).^2).*exp(-1e3*(t-.2).^2);

% plot
figure;
plot(t,g);
xlabel('t [s]');ylabel('g(t)');
title('true sparse signal')

figure;
plot(t,w);
xlabel('t [s]');ylabel('w(t)');
title('band pass filter')


% fourier transform of w
wf = fft(w);

% SPOT operator to perform convolution.
C = opDFT(N)'*opDiag(wf)*opDFT(N);
f = C*g;

% plot
figure;
plot(t,f);
xlabel('t [s]');ylabel('f(t)');
title('convolution result')




%% spgl1 and pqnl1
cond(full(C))
    %% lasso
    opts.iterations = 100;
    tau = norm(g,1);
      
    opts.fid = fopen('spg','w');
    [x_spg,r_spg,g_spg,info_spg] = spgl1(C, f, tau, [], zeros(size(g)), opts);
    
    opts.fid = fopen('pqn','w');
    [x_pqn,r_pqn,g_pqn,info_pqn] = pqnl1_2(C, f, tau, [], zeros(size(g)), opts);

    h = figure; 
    subplot(3,1,1); plot(g); title('original sparse signal')
    subplot(3,1,2); plot(x_spg);title('x_spg')
    subplot(3,1,3); plot(x_pqn);title('x_pqn')
    save(h,'deconvolution result')
  



    

