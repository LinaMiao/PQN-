function [f ky kh h] =  kykh(t,xr,xs,squeeze)
% function [f ky kh h] =  kykh(t,xr,xs)
% input :  t - time
%       :  xr- receiver
%       :  xs- source
% output : f - frequency
%        : ky - midpoint wavenumber
%        : kh - offset wavenumber
%        : h  - offset coordinate


% only works if even length of xr


if not(exist('squeeze','var'))
    squeeze = 0;
end

n = length(xr);
[xxr xxs] = ndgrid(xr,xs);
hh = .5*(xxs - xxr);
yy = .5*(xxs + xxr);


% offset coordinate
h = unique(hh)';
dh = h(2)-h(1);
h = [h,h(end)+dh];

% midpoint coordinate
y = unique(yy)';
dy = y(2)-y(1);
y = [y,y(end)+dy];

% squeeze
if squeeze == 1
    y = y(2:2:end);
end

% compute corresponding frequency and wave number
[f, ky ,kh] = fkk(t,y,h);

