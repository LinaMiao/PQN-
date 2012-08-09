function data = DSR_mig_time_inv(image,t,x,y,z,v,option)
% inverse dsr extrapolation if input in time domain
% function output = test11(data,t,xr,xs,z,v,option)
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
if option == 0 % migration
    image = reshape(image,length(z),length(x));
    eximage = zeros(length(z),length(x),length(y));
    % adjoint of taking diag if migration
    for iz = 1:length(z)
        D = diag(ones(length(image(iz,:)),1));
        eximage(iz,:) = opDiag(vec(D))'*vec(diag(vec(image(iz,:))));
    end
    clear image;
    image = eximage;
else
    image = reshape(image,length(z),length(x),length(y));
end

    

% initialize output
data = zeros(length(t),length(x),length(y));

% depth step 
dz = z(2) - z(1);

% f-k-k grid
[f,kx,ky] = fkk(t,x,y);    
[ff,kkx,kky] = ndgrid(f,kx,ky);

% the summing operator
S = opKron(opDirac(length(ky)),opDirac(length(kx)),ones(1,length(f)));

% loop over depth levels
for iz = length(z):-1:1
    % adjoint of summing over all frequency
    ads = S'*vec(image(iz,:));
    if iz == length(z)
        temp = ads;
    else
        temp = fft(data,[],1)/sqrt(size(data,1)) + reshape(ads,size(data));
    end
    % adjoint of survey sinking
    data = DSR_inv_step(temp,ff,kkx,kky,dz,v(iz));
end






function v = DSR_inv_step(u,ff,kkx,kky,dz,v)

% inverse DSR 
% input f-x-x
% output t-x-x

if size(size(u)) ~= 3
    u = reshape(u,size(ff));
end

Px = 2*pi*sqrt((ff/v).^2-kkx.^2);
Py = 2*pi*sqrt((ff/v).^2-kky.^2);

Px = real(Px)+1i*abs(imag(Px));
Py = real(Py)+1i*abs(imag(Py));

E = exp(1i*abs(dz)*(Px + Py));


v = u;
v = fft(v,[],2)/sqrt(size(u,2));
v = fft(v,[],3)/sqrt(size(u,3));

% apply the adjoint of phase shift
E =  exp(1i*abs(dz)*(Px + Py));
v = opDiag(vec(E))'*vec(v);

% reshape after spot operation
v = reshape(v,size(ff));

v = ifft(v,[],1)*sqrt(size(v,1));
v = ifft(v,[],2)*sqrt(size(v,2));
v = ifft(v,[],3)*sqrt(size(v,3));









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





