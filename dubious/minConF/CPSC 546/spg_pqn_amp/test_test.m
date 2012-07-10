clear;

delta = 1/2;
rho   = 1/20;
 
nscale  = 1;
scale0  = 8;

for iscale=1:nscale
    N = 10^3*(scale0+iscale);
    n = floor(delta*N); 
    k = floor(rho*n);
    inds=randperm(N);
    ind=inds(1:n);
    R=opRestriction(N,ind);
    A=1/sqrt(n)*R*sqrt(N)*opRomberg([N]);
    %A = double(A);
    %A = rand(n,N); 
    x = [sign(rand(k,1) - 0.5); zeros(N-k,1)];
    x0 = [sign(rand(k,1) - 0.5); zeros(N-k,1)]; 
    x = x(randperm(N));
    tau0=norm(x,1);
    b = A*x;
    opts = spgSetParms('optTol',1e-3,'verbosity',1,'iterations',500);
    [xspg,rspg,gspg,infospg{iscale}] = spgl1(A, b, 0, 1e-4, zeros(size(A,2),1), opts);
    [xpqn,rpqn,gpqn,infopqn{iscale}] = pqnl1(A, b, 0, 1e-4, zeros(size(A,2),1), opts);
    tampstart = toc;
    [xamp,infoamp{iscale}] = reconstructAmpcountfast(A, b, infospg{iscale}.iter, 1e-2);
     
    timeamp(iscale) = toc-tampstart;
    timespg(iscale) = infospg{iscale}.timeTotal;
    timepqn(iscale) = infopqn{iscale}.timeTotal;
    
end
save test xspg xpqn xamp N nscale infospg infopqn infoamp timespg timepqn timeamp
%load test;

figure;plot(1:nscale,timespg,1:nscale,timepqn,1:nscale,timeamp);legend('spg','pqn','amp')

%%
figure('Name','Solution paths')
plot(infospg{end,end}.xNorm1,infospg{end,end}.rNorm2,infopqn{end,end}.xNorm1,infopqn{end,end}.rNorm2,infoamp{end,end}.xNorm1,infoamp{end,end}.rNorm2);hold on
hold on;
scatter(infospg{end,end}.xNorm1,infospg{end,end}.rNorm2);
scatter(infopqn{end,end}.xNorm1,infopqn{end,end}.rNorm2);
scatter(infoamp{end,end}.xNorm1,infoamp{end,end}.rNorm2);hold off
xlabel('one-norm model')
ylabel('two-norm residual')
title('Solutions paths')
legend('SPGl1','PQNl1','AMP')
axis tight