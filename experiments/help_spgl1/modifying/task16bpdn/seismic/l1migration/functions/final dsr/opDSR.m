function op = opDSR(dim,t,xr,xs,z,v,option)

% by default, extended migration for input in frequency domain
if not(exist('option','var'))
    option.freq = 1;
    fprintf('not specify option, by default extended migration for input in frequency domain')
end

% Determine sizes
m = prod(dim);

[xxr xxs] = ndgrid(xr,xs);
hh = .5*(xxs-xxr);
h = unique(hh)';
[f,kh] = fkk(t,h);


[ff,kkh] = ndgrid(f,kh);
p = kkh./ff;
p = unique(p);

n = length(xs)*length(z)*length(p);

% Construct the operator
fh = @(data,mode) opDSRmig_intrnl(data, t, xr, xs, z, v ,option,mode);
op = opFunction(n, m, fh);


function y = opDSRmig_intrnl(data,t,xr,xs,z,v,option,mode)

    if (mode == 0)
        
        y = {m,m,[0,0,0,0],{'opDSR'}};
        
    elseif (mode == 1)
        
        y = DSR(data,t,xr,xs,z,v,option);
        y = vec(y);
        
    else % mode = -1
        
        y = inv_DSR(data,t,xr,xs,z,v,option);
        y= vec(y);
        
    end % end of if mode == 0




end

end