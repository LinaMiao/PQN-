function [f kx ky] = fkk(t,x,y)
% compute f kx ky via t x y
if (exist('y','var'))
    dy   = y(2) - y(1);
    ny   = length(y);
    fny  = .5/dy;
    if mod(ny,2) == 1
        df   = 2*fny/(ny-1);
        ky   = [0:df:fny -fny:df:-df];
    else
        df   = 2*fny/ny;
        ky   = [0:df:fny -fny+df:df:-df];
    end
end

if (exist('x','var'))
    dx   = x(2) - x(1);
    nx   = length(x);
    fny  = .5/dx;
    if mod(nx,2) == 1
        df   = 2*fny/(nx-1);
        kx   = [0:df:fny -fny:df:-df];
    else
        df   = 2*fny/nx;
        kx   = [0:df:fny -fny+df:df:-df];
    end
end

if (exist('t','var'))
    dt   = t(2)-t(1);
    nt   = length(t);
    fny  = .5/dt;
    if mod(nt,2) == 1
        df   = 2*fny/(nt-1);
        f    = [0:df:fny -fny:df:-df];
    else
        df   = 2*fny/nt;
        f    = [0:df:fny -fny+df:df:-df];
    end
end

    
    
    

  




