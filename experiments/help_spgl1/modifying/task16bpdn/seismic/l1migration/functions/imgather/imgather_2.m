function [output] = imgather_2(input,z,xr,xs,t)

% function [output] = imgather_2(input,z,xr,xs,t)
% input :  input(m,z,p); 
% output :  output(z,xr,xs);

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
if size(size(input)) ~= 3
    input = reshape(input,nm,nz,np);
end

% initialize output
output = zeros(nz,nr,nm);


% for each midpoint
for i = 1:nm 
    
    % adjoint of radon for this particular midpoint gather
    d = squeeze(input(i,:,:));
    hi = hh(find(mm == m(i)))'; 
    temp = lpradon(d,t,hi,p,1,-1);
    
    % for each depth level put output of adjoint radon to right position
    outputi = zeros(nz,nr*nm);
    for j = 1:nz
        outputi(j,find(mm == m(i))) = temp(j,:);
    end
    tt = reshape(outputi,nz,nr,nm);
    output = output + tt;
end

    

    


