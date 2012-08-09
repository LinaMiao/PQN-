function output = lpradon(input,t,h,q,power,mode,fmin,fmax)

% fmin and fmax is user defined freq range in Hz

% if fmin and fmax is not given, determine here
fny = .5/(t(2)-t(1));
if not(exist('fmin','var'))
    fmin = min(2,fny);
end

if not(exist('fmax','var'))
    fmax = min(fny,25);
end

if fmin == fmax
    fmin = 1/t(end);
end

% transform fmin and fmax to 'circular frequency'
fmin = 2*pi*fmin;
fmax = 2*pi*fmax;


if mode == -1
    output = inverse_radon_freq(input,t,h,q,power,fmin,fmax);
else
    output = forward_radon_freq(input,t,h,q,power,fmin,fmax);
end

end

function [m] = forward_radon_freq(d,t,h,q,N,fmin,fmax)
% input d(t,h), time and offset
% output m(t,p), time and waveparameter

% determine sizes
[nt,nh] = size(d);
nq = max(size(q));

% forward fourier on the first dimension of d->D(t,h)
D = fft(d,[],1)/sqrt(size(d,1));

% compute frequency coordinate
faxis = 2.*pi*fkk(t);

% initialize output
M = zeros(nt,nq);

% imaginary unit
i = sqrt(-1);

% window the input, only radon on a cirtain frequency range
% plausible frequency range
df = 2*pi*(1/nt);
idxmin = find((faxis>fmin-df/2).*(faxis<=fmin+df/2),1);
idxmax = find((faxis>fmax-df/2).*(faxis<=fmax+df/2),1);
% mask matrix
mask = zeros(size(D));
mask(idxmin:idxmax,:) = 1;
% apply mask
D = mask.*D;

% Radon, for each f and p, sum over all the h   
for ifreq=1:size(D,1)
    f = faxis(ifreq);
    L = exp(i*f*(h.^N)'*q);
    y = D(ifreq,:)';
    xa = L'*y;
    x  =  xa;
    M(ifreq,:) = x';
end

% inverse fourier transform on the first dimension M(f,p)->m(t,p)
m = (ifft(M,[],1))*sqrt(size(M,1));

end

function [d]=inverse_radon_freq(m,t,h,p,N,fmin,fmax)
% input m(t,p), time and waveparameter
% output d(t,h), time and offset

% determine sizes
[nt,nq] = size(m);
nh = length(h);

% for N=2 Sacchi is doing offset expension, which I did not keep, if
% necessary, check in SeisLab code

% forward fourier on the first dimension of d->D(k,h)
M = fft(m,[],1)/sqrt(size(m,1));

% compute frequency coordinate
faxis = 2.*pi*fkk(t);

% initialize output
D = zeros(nt,nh);

% imaginary unit
i = sqrt(-1);

% apply the adjoint of mask used in forward function
df = 2*pi*(1/nt);
idxmin = find((faxis>fmin-df/2).*(faxis<=fmin+df/2),1);
idxmax = find((faxis>fmax-df/2).*(faxis<=fmax+df/2),1);
% mask matrix
mask = zeros(size(M));
mask(idxmin:idxmax,:) = 1;
% apply adjoint mask
M = opDiag(mask)'*vec(M);
M = reshape(M,size(m));

% Radon, for each f and h, sum over all the p   
for ifreq=1:size(M,1)
    f = faxis(ifreq);
    L = exp(i*f*(h.^N)'*p);
    x = M(ifreq,:)';
    y = L * x;
    D(ifreq,:) = y';
end

% inverse fourier transform on the first dimension M(f,p)->m(t,p)
d = (ifft(D,[],1))*sqrt(size(D,1));

end