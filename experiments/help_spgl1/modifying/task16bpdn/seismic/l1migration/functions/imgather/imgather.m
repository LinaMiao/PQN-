function [output] = imgather(data,z,xr,xs,t)

% function [output] = imgather(data,z,xr,xs,t)
% input :  data(z,xr,xs); 
% output :  output(m,z,p);

% compute all the offset coordinate
[xxr xxs] = ndgrid(xr,xs);
mm = .5*(xxs+xxr);
hh = .5*(xxs-xxr);
% m = unique(mm)';
h = unique(hh)';
[f,kh] = fkk(t,h);

% take midpoints corresponding to source location
m = xs; 

% compute waveparameter p 
[ff,kkh] = ndgrid(f,kh);
p = kkh./ff;
pmin = min(min(isfinite(p))); pmax = max(max(isfinite(p)));p = unique(p);
np = length(p);
p = linspace(pmin,pmax,np);

%determine sizes
nm = length(m);
nh = length(h);
nz = length(z);
nr = length(xr);
ns = length(xs);
nt = length(t);

% reshape input
if size(size(data)) ~= 3
    data = reshape(data,nz,nr,ns);
end


% initialize output
output = zeros(round(nm),nz,np);

% for each midpoint
for i = 1:nm 
    hi = hh(find(mm == m(i)))'; 
    d = zeros(nz,length(hi)); % form a common mid point gather d
    for j = 1:nz
        temp = squeeze(data(j,:,:));
        d(j,:) = temp(find(mm == m(i)));
    end
       
    % apply radon to each common offset gather
    output(i,:,:) = lpradon(d,t,hi,p,1,1);
end


   
    


