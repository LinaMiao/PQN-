function output = DSR_mig_time(data,t,x,y,z,v,option)
% dsr extrapolation if input in time domain
% function output = DSR_time(data,t,xr,xs,z,v,option)
% x is xr -- receiver grid
% y is xs -- source grid
% t       -- time axis
% z       -- vertical axis
% v       -- velocity as a function of z, must be the same size of z
% option  -- 0:migration
%         -- 1:extended migration


if norm(x-y)
    warning('source and receiver grids must be the same!');
end

if length(v) - length(z)
    warning('velocity must be the same size as z')
end

% reshape input
if size(size(data)) ~= 3
    data = reshape(data,length(data)/(length(x)*length(y)),length(x),length(y));
end

% initialize image
if option == 0
    output = zeros(length(z),length(x));
else 
    output = zeros(length(z),length(x),length(x));
end

% depth step 
dz = z(2) - z(1);

% basic computation of frequencies and wave numbers
[f kx ky] =  fkk(t,x,y);

% f-k-k grid
[ff,kkx,kky] = ndgrid(f,kx,ky);

% loop over depth levels
for iz = 1:length(z)
    E = DSR_step(data,ff,kkx,kky,dz,v(iz));
    % E = reshape(E,size(data));
    data = ifft(E,[],1)*sqrt(size(E,1));
    if option == 0 % migration
        output(iz,:) = diag(squeeze(sum(E,1)));
    else
        output(iz,:,:) = squeeze(sum(E,1));
    end
end




function v = DSR_step(u,ff,kkx,kky,dz,v)

% DSR 
% input t-x-x
% output f-x-x

if size(size(u)) ~= 3
    u = reshape(u,size(ff));
end

Px = 2*pi*sqrt((ff/v).^2-kkx.^2);
Py = 2*pi*sqrt((ff/v).^2-kky.^2);

Px = real(Px)+1i*abs(imag(Px));
Py = real(Py)+1i*abs(imag(Py));

E = exp(1i*abs(dz)*(Px + Py));


v = fft(u,[],1)/sqrt(size(u,1));
v = fft(v,[],2)/sqrt(size(u,2));
v = fft(v,[],3)/sqrt(size(u,3));

v = E.*v;

v = ifft(v,[],2)*sqrt(size(v,2));
v = ifft(v,[],3)*sqrt(size(v,3));


% F1 = opFFT3(size(ff));
% F2 = opFFT2(size(ff));
% v = F2'*opDiag(vec(E))*F1*vec(u);






function [f kx ky] = fkk(t,x,y)
% compute f kx ky via t x y
dt   = t(2)-t(1);
nt   = length(t);
dx   = x(2) - x(1);
nx   = length(x);
dy   = y(2) - y(1);
ny   = length(y);
fny  = .5/dt;
if mod(nt,2) == 1
    df   = 2*fny/(nt-1);
    f    = [0:df:fny -fny:df:-df];
else
    df   = 2*fny/nt;
    f    = [0:df:fny -fny+df:df:-df];
end
fny  = .5/dx;
if mod(nx,2) == 1
    df   = 2*fny/(nx-1);
    kx   = [0:df:fny -fny:df:-df];
else
    df   = 2*fny/nx;
    kx   = [0:df:fny -fny+df:df:-df];
end
fny  = .5/dy;
if mod(ny,2) == 1
    df   = 2*fny/(ny-1);
    ky   = [0:df:fny -fny:df:-df];
else
    df   = 2*fny/ny;
    ky   = [0:df:fny -fny+df:df:-df];
end





