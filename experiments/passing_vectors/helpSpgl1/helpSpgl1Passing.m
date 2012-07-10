% this experiment is to test whether pqnl1 can work for the example given
% by help spgl1

% pqnl1 is used as passing vectors version, whether passing works and
% whether it improves performance


%% problem setting
% m = 240; n = 1024; k = 40; % m rows, n cols, k nonzeros.
% p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
% A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
% b  = A*x0 + 0.005 * randn(m,1);
% 
% opts.optTol = 1e-4;
% opts.iterations = 200;
% 
% save passing m n k p A b opts


%%
clear;
load passing
%% spgl1
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, b, 0, 1e-3, [], opts); % Find BP sol'n.

%% pqnl1_2
[x_pqn,r_pqn,g_pqn,info_pqn] = pqnl1_2(A, b, 0, 1e-3, zeros(size(A,2),1), opts); % Find BP sol'n.

%% pqnl1_passing
[x_passing,r_passing,g_passing,info_passing] = pqnl1_passing(A, b, 0, 1e-3, zeros(size(A,2),1), opts); % Find BP sol'n.

%% show result
figure; subplot(3,1,1);plot(x_spg);subplot(3,1,2);plot(x_pqn);subplot(3,1,3);plot(x_passing);

figure('Name','Solution paths')
plot(info_spg.xNorm1,info_spg.rNorm2,info_pqn.xNorm1,info_pqn.rNorm2,info_passing.xNorm1,info_passing.rNorm2);hold on
scatter(info_spg.xNorm1,info_spg.rNorm2);
scatter(info_pqn.xNorm1,info_pqn.rNorm2);
scatter(info_passing.xNorm1,info_passing.rNorm2);hold off
xlabel('one-norm model')
ylabel('two-norm residual')
title('Solutions paths')
legend('SPGl1','PQNl1','PQN_passing')
axis tight

info_spg
info_pqn
info_passing

%% check functions
open ./minConF_PQN_passing.m
open ./pqnl1_passing.m
open ./minConF_PQN_2.m
open ./pqnl1_2.m