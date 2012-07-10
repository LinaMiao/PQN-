% This is for the second problem in homework2 from CPSC 546
% use Newton's method to find the reciprocal for a number d

function [reciprocal,itn] = newreci(d,x0)
f = @(x) log(x) + log(d); 
x = x0;
tol = 10^-6;
itn = 0;
while abs(f(x))>=tol
    x = x-(log(x)+log(d))*x
    itn = itn+1
end
reciprocal = x;

% with differents initial guess, Newton's method converges with different
% rates.
