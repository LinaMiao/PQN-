function CIG = DSR(data,t,xr,xs,z,v,option)

% form extended image gather
if option.freq == 1
    f = fkk(t);
    DSR_step = opDSR_mig_freq([length(f),length(xr),length(xs)],f,xr,xs,z,v(1:length(z)),1);
else
    DSR_step = opDSR_mig_time([length(t),length(xr),length(xs)],t,xr,xs,z,v(1:length(z)),1);
end
% image = DSR_step*vec(data);

% form common waveparameter gather via Radon
dim =  [length(z),length(xr),length(xs)];
op = opImgather(dim,z,xr,xs,t);
CIG = op*DSR_step*vec(data);

