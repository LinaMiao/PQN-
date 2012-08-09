function [output] = mixNorm(X,p,q)
% compute the mixed norm of a matrix
% input : X p q 
%       || X ||_p,q = ( sum(j:1-n) || X^j->||^p_q ) ^1/p
%       X^j-> : is the vector whose entries form the jth row of X


if nargin < 3
    error('not enough input')
end

n = size(X,1);
output = 0;
for j = 1: n
    xj = vec(X(j,:));
    output = output + norm(xj,q)^p;
end
output = output^(1/p);


    
    