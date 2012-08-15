%% addpath for PQN working
cd ../../../../functions;
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task12illconditioned
%stream = RandStream.getGlobalStream;
%reset(stream);


%% sample matrix and options
i = 1;
for m = [ 200 250 300 350 400 450 500]
    %m = 200; 
    n = 512; k = 20; % m rows, n cols, k nonzeros.
    A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
    [u s v] = svd(A); 

    %figure;plot(diag(s));title('singular values of A')
    ns = length(diag(s));
    nn = linspace(0,1,ns);
    s_ill = exp(-nn.^.1);
    s_ill = s_ill - (1-1e-6)*min(s_ill);
    condition_number = max(s_ill)/min(s_ill)

    s_new = zeros(m,n);
    s_new(1:min(m,n),1:min(m,n)) = diag(s_ill);


    %figure;plot(diag(s_new));title('proposed singular values')
    A_ill = u'*s_new*v;

    opts.iterations = 100;
    %opts.verbosity = 0;

    % save temp A m n k opts
    % clear;
    % load temp

    %% problem setting
    p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
    %figure;plot(x0)
    b  = A_ill*x0;

    tau = norm(x0,1);

    %% Lasso
    opts.fid = fopen(strcat(['spg' num2str(i)]),'w');
    [x_spg,r_spg,g_spg,info_spg] = spgl1(A_ill, b, tau, [], zeros(size(x0)), opts); 
    
    opts.fid = fopen(strcat(['pqn ' num2str(i)]),'w');
    [x_pqn,r_pqn,g_pqn,info_pqn] = pqnl1_2(A_ill, b, tau, [], zeros(size(x0)), opts);

    h = figure(i);
    subplot(2,1,1); plot(x_spg);title('x_spg')
    subplot(2,1,2); plot(x_pqn);title('x_pqn')
    saveas(h,strcat([' ' num2str(i)]));

    i = i + 1;

end