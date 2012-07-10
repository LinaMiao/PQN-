function [ v ] = randvector( n, meanv, stdv, norm_dist)

 

    

%RADMVECTOR Summary of this function goes here

%   Detailed explanation goes here

    if nargin < 4

        norm_dist = 1;

    end

    

    if norm_dist

        v = randn(n, 1);

    else

        v = rand(n, 1);

    end

    

    

    v = v + meanv - mean(v);

    v = (v-meanv) / std(v) * stdv + meanv;

    

    

 

end

 



