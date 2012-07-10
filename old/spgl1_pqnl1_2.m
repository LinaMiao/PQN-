load test
%opts = spgSetParms('optTol',1e-4,'iterations',200);
opts.optTol = 1e-4;
opts.iterations = 200;
%% spgl1
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, b, 0, 1e-3, [], opts); % Find BP sol'n.

%% pqnl1
[x_pqn,r_pqn,g_pqn,info_pqn] = pqnl1_2(A, b, 0, 1e-3, zeros(size(A,2),1), opts); % Find BP sol'n.
