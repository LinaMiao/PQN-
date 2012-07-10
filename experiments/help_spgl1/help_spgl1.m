% this experiment is to test whether pqnl1 can work for the example given
% by help spgl1

%% problem setting
m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
b  = A*x0 + 0.005 * randn(m,1);

opts.optTol = 1e-4;
opts.iterations = 200;
%% spgl1
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, b, 0, 1e-3, [], opts); % Find BP sol'n.

%% pqnl1
[x_pqn,r_pqn,g_pqn,info_pqn] = pqnl1_2(A, b, 0, 1e-3, zeros(size(A,2),1), opts); % Find BP sol'n.

%% show result
figure; subplot(2,1,1);plot(x_spg);subplot(2,1,2);plot(x_pqn);
info_spg
info_pqn

%% check functions
open ./minConF_PQN_2.m
open ./pqnl1_2.m