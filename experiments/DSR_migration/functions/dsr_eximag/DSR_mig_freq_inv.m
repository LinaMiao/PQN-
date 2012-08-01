function data = DSR_mig_freq_inv(image,f,x,y,z,v,option)

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
data = zeros(length(f),length(x),length(y));

% depth step 
dz = z(2) - z(1);

% f-k-k grid
[kx,ky] = fkk(x,y);    
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
        temp = data + reshape(ads,size(data));
    end
    % adjoint of survey sinking
    data = DSR_inv_step(temp,ff,kkx,kky,dz,v(iz));
end



function v = DSR_inv_step(u,ff,kkx,kky,dz,v)
% the adjoint of DSR_step
% input : u(f,x,x),ff,kkx,kky,dz,v
% output : v(f,x,x) 

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

% apply the adjoint of phase shift
E =  exp(1i*abs(dz)*(Px + Py));
spec = opDiag(vec(E))'*vec(spec);

% reshape after spot operation
spec = reshape(spec,size(ff));

% inverse fkk transform
v = spec;
v = ifft(v,[],2)*sqrt(size(v,2));
v = ifft(v,[],3)*sqrt(size(v,3));


