function image = DSR_mig_freq(data,f,x,y,z,v,option)

% DSR migration if data is given in f-x-x domain

% x is xr -- receiver grid
% y is xs -- source grid
% f       -- frequencies axis
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
    image = zeros(length(z),length(x));
else 
    image = zeros(length(z),length(x),length(x));
end

% depth step 
dz = z(2) - z(1);

% f-k-k grid
[~,kx,ky] = fkk(f,x,y); 
if size(data,1) ~= length(f) % only positive frequency is provided
    n = length(f);
    if mod(n,2) == 1
       tmp = [f , -(f(end:-1:2))];
    else
       tmp = [f , -(f(end-1:-1:2))];
    end
    f = tmp;
end
    
[ff,kkx,kky] = ndgrid(f,kx,ky);

% loop over depth levels
for iz = 1:length(z)
    % survey sinking
    E = DSR_step(data,ff,kkx,kky,dz,v(iz));
    % updata data for each depth level
    data = E;
    % imaging condition
    if option == 0
        image(iz,:) = diag(squeeze(sum(E,1)));
    else
        image(iz,:,:) = (squeeze(sum(E,1)));
    end
end


function v = DSR_step(u,ff,kkx,kky,dz,v)

% reshape input
if size(size(u)) ~= 3
    u = reshape(u,size(ff));
end

% f-k-k transform
spec = u;
spec = fft(spec,[],2)/sqrt(size(spec,2));
spec = fft(spec,[],3)/sqrt(size(spec,3));

% DSR operators
Px = 2*pi*sqrt((ff/v).^2-kkx.^2);
Py = 2*pi*sqrt((ff/v).^2-kky.^2);

Px = real(Px)+1i*abs(imag(Px));
Py = real(Py)+1i*abs(imag(Py));

% apply phase shift
spec = exp(1i*abs(dz)*(Px + Py)).*spec;

% inverse fkk transform
v = spec;
v = ifft(v,[],2)*sqrt(size(v,2));
v = ifft(v,[],3)*sqrt(size(v,3));

% take real part
% v = real(v);

