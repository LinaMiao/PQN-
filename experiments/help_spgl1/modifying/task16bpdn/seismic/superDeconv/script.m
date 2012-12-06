close all
clear
rand('state',271)

% spike signal on 0 to 1
n = 1024;
dt = 1/(n-1);
t = 0:dt:1;
 

fc = 128; 
dx = ceil(2/fc * n); 
kmax = floor(n/dx);
k = min(floor(fc/5),kmax); 
temp = randperm(kmax);
 
p = temp(1:k) .* dx;
dynamic_range = 0;
x0 = zeros(n,1); 
x0(p) = sign(randn(k,1)) .* (1 + 10.^(rand(k,1).*(dynamic_range/20))); % amplitudes can also be complex

max(abs(x0(p)))/min(abs(x0(p)))
h = figure; stem(t,x0); title('original signal')
axis('tight'); saveas(h,'original_signal')
 
 
% low pass fourier filter in freq domain and in time domain
f = zeros(n,1);
f(1) = 1;
ff = fft(f,n);
ff(fc : n-fc) = 0;
f_low = ifft(ff);

 
%time-domain equivalent (convolution operator formation)
[M N]=size(x0);
B=fc;
ts=t-0.1*max(t);
Kf=2.*B.*sinc(2.*B.*ts); %sinc function
S=opConvolve(M,N,Kf');
test3=S*x0;
 
% % test of operator
% test1 = ifft(fft(x0)/sqrt(n) .* (ff/sqrt(n)))*sqrt(n);
% figure;
% subplot(2,1,1);plot(real(test1))
% subplot(2,1,2);plot(real(test3))
 
h = figure;
plot(Kf);title('source_wavelet')
saveas(h,'source_wavelet')
 
 
% solve BP problem with spgl1
options.iterations=103;
options.fid = fopen('spg.txt','w');
A=S;
b=test3;
[x_td,r,g,info_spg] = spgl1(A, b, 0, 0, [],options);%time domain reconsruction
info_spg.rNorm
save infospg info_spg
h = figure;
subplot(2,1,1);stem(x0);axis tight; ylim([min(x0) max(x0)])
subplot(2,1,2);stem(x_td); axis tight;ylim([min(x0) max(x0)])
title('spg')
saveas(h,'spg');

sigma_ref = info_spg.rNorm;
options.fid = fopen('pqn.txt','w');

[x_td_2,r,g,info_pqn] = pqnl1_2(A, b, 0, 0,zeros(length(x_td),1),options,sigma_ref);%time domain reconsruction
h = figure;
subplot(2,1,1);stem(x0);axis tight; ylim([min(x0) max(x0)])
subplot(2,1,2);stem(x_td_2); axis tight; ylim([min(x0) max(x0)])
title('pqn')
saveas(h,'pqn')
save infopqn info_pqn


h = figure('Name','Solution paths');
plot(info_spg.xNorm1,info_spg.rNorm2,info_pqn.xNorm1,info_pqn.rNorm2);hold on
scatter(info_spg.xNorm1,info_spg.rNorm2);
scatter(info_pqn.xNorm1,info_pqn.rNorm2);hold off
legend('spg','pqn')
axis tight
saveas(h,'solution_path');

