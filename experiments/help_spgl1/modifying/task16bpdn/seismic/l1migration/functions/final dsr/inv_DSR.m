function data = inv_DSR(CIG,t,xr,xs,z,v,option)

% form extended image
dim =  [length(z),length(xr),length(xs)];
op = opImgather(dim,z,xr,xs,t);
% Image = op'*vec(CIG);

% back propagate via inv DSR
if option.freq == 1
    f = fkk(t);
    DSR_step = opDSR_mig_freq([length(f),length(xr),length(xs)],f,xr,xs,z,v(1:length(z)),1);
else
    DSR_step = opDSR_mig_time([length(t),length(xr),length(xs)],t,xr,xs,z,v(1:length(z)),1);
end

data = DSR_step'*op'*vec(CIG);





