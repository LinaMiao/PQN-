%% add working path of tools and functions
%addpath(genpath('e:\research\Tools\spot-slim'))
% addpath(genpath('e:\research\Tools\pSPOT'))
cd ../../../../../pqnl1;
addpath(genpath(pwd))
rmpath(genpath('./minConF'))
cd ../experiments/help_spgl1/modifying/task12illconditioned
addpath(genpath(pwd))
cd ./convolution
rmpath('/Volumes/Users/linamiao/Dropbox/PQN/pqnl1/minConF/')


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
    %% lasso
    opts.iterations = 100;
    tau = norm(g,1);

    [x_spg,r_spg,g_spg,info_spg] = spgl1(C, f, tau, [], zeros(size(g)), opts);
    opts.iterations = 50;
    [x_pqn,r_pqn,g_pqn,info_pqn] = pqnl1_2(C, f, tau, [], zeros(size(g)), opts);

    figure; 
    subplot(3,1,1); plot(g); title('original sparse signal')
    subplot(3,1,2); plot(x_spg);title('x_spg')
    subplot(3,1,3); plot(x_pqn);title('x_pqn')
    
    %% BPDN
    % noisy signal
    f = C*g+ 1e-2*randn(N,1);

    % plot
    figure;
    plot(t,f);
    xlabel('t [s]');ylabel('f(t)');
    title('convolution result')
    
    opts.iterations = 100;
    

    [x_spg,r_spg,g_spg,info_spg] = spgl1(C, f, 0, 1e-2, zeros(size(g)), opts);
    [x_pqn,r_pqn,g_pqn,info_pqn] = pqnl1_2(C, f, 0, 1-2, zeros(size(g)), opts);

    figure; 
    subplot(3,1,1); plot(g); title('original sparse signal')
    subplot(3,1,2); plot(x_spg);title('x_spg')
    subplot(3,1,3); plot(x_pqn);title('x_pqn')

    figure('Name','Solution paths')
    plot(info_spg.xNorm1,info_spg.rNorm2,info_pqn.xNorm1,info_pqn.rNorm2);hold on
    scatter(info_spg.xNorm1,info_spg.rNorm2);
    scatter(info_pqn.xNorm1,info_pqn.rNorm2);hold off
    legend('spg','pqn')
    axis tight



    

