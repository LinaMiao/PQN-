% this script is to test whether current version of pqnl1 can solve problem
% in spgl1's example.

% addpath for pqnl1_2 working
addpath(genpath('./pqnl1')); 

%% problem setting
m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
b  = A*x0 + 0.005 * randn(m,1);

save test m n k p x0 A b;