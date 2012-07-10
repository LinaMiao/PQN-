% In this script we investigate the timing as a function of the numer of
% matvecs as the system size increases
clear;%close all;

delta = 1/200;
rho   = 1/10;

nscale  = 2;
scale0  = 11;
tol     = 1e-3;

nmatvec = 1;
matrixvector0    = 200;
matrixvector1    = 200;
matrixvectors     = floor(linspace(matrixvector0,matrixvector1,nmatvec));

errorspg=zeros(nmatvec,nscale,5);
erroramp=zeros(nmatvec,nscale,5);

if 1,
    for iscale=1:nscale
        N = 2^(scale0+iscale);
        n = floor(delta*N);
        k = floor(rho*n);
        inds=randperm(N);
        ind=inds(1:n);
        R=opRestriction(N,ind);
        A=1/sqrt(n)*R*sqrt(N)*opRomberg([N]);
        x = [sign(rand(k,1) - 0.5); zeros(N-k,1)];
        x = x(randperm(N));
        tau0=norm(x,1);
        % Generate Measurements
        b = A*x;
        for imatvec = 1:nmatvec
            % Do the one-norm solve
            nmatrixvectors = matrixvectors(imatvec);
            opt.verbosity = 0;
            opt.maxMatvec = nmatrixvectors;
            %[xhatspg,r,~,infospg{irea,iscale}] = spg_bpdn(A, b, tol*norm(b,2),opt);
            [xhatspg,r,~,infospg{imatvec,iscale}] = localspgl1(A, b, 0, 1e-2, zeros(size(A,2),1), opt);
            errorspg(imatvec,iscale,1)=N;
            errorspg(imatvec,iscale,2)=infospg{imatvec,iscale}.nProdA+infospg{imatvec,iscale}.nProdAt;
            errorspg(imatvec,iscale,3)=norm(r,2);
            errorspg(imatvec,iscale,4)=norm(x-xhatspg,2);
            errorspg(imatvec,iscale,5)=norm(x,2);
            errorspg(imatvec,iscale,6)=infospg{imatvec,iscale}.iter;
            
            
            [xhatpqn,r,~,infopqn{imatvec,iscale}] = pqnl1(A, b, 0, 1e-1, zeros(size(A,2),1), opt);
            errorpqn(imatvec,iscale,1)=N;
            errorpqn(imatvec,iscale,2)=infopqn{imatvec,iscale}.nProdA+infopqn{imatvec,iscale}.nProdAt;
            errorpqn(imatvec,iscale,3)=norm(r,2);
            errorpqn(imatvec,iscale,4)=norm(x-xhatpqn,2);
            errorpqn(imatvec,iscale,5)=norm(x,2);
            errorpqn(imatvec,iscale,6)=infopqn{imatvec,iscale}.iter;
           
%             [xhatamp,infoamp{imatvec,iscale}] = reconstructAmpcountfast(A, b, niterations, tol);
%             erroramp(imatvec,iscale,1)=N;
%             erroramp(imatvec,iscale,2)=infoamp{imatvec,iscale}.nProdAAt;
%             erroramp(imatvec,iscale,3)=norm(A*xhatamp-b,2);
%             erroramp(imatvec,iscale,4)=norm(x-xhatamp,2);
%             erroramp(imatvec,iscale,5)=norm(x,2);
%             erroramp(imatvec,iscale,6)=infoamp{imatvec,iscale}.iter;
            
        end
    end
    
    save varyingiterations errorspg errorpqn erroramp N delta rho matrixvectors infospg infopqn %infoamp
else
    load varyingiterations
end
%

%%
figure;plot(squeeze(errorspg(1,:,1)),squeeze(errorspg(:,:,4))./squeeze(errorspg(:,:,5)),'-','LineWidth',2);title('SPGl1')
xlabel('N');ylabel('relative error');% legend(num2str(matrixvectors(1)),num2str(matrixvectors(2)),num2str(matrixvectors(3)))%,num2str(iterations(4)),num2str(iterations(5)),num2str(iterations(6)),num2str(iterations(7)),num2str(iterations(8)),num2str(iterations(9)),num2str(iterations(10)))
axis tight
figure;plot(squeeze(errorpqn(1,:,1)),squeeze(errorpqn(:,:,4))./squeeze(errorpqn(:,:,5)),'-','LineWidth',2);title('pqnl1')
xlabel('N');ylabel('relative error');% legend(num2str(matrixvectors(1)),num2str(matrixvectors(2)),num2str(matrixvectors(3)))%,num2str(iterations(4)),num2str(iterations(5)),num2str(iterations(6)),num2str(iterations(7)),num2str(iterations(8)),num2str(iterations(9)),num2str(iterations(10)))
axis tight
%%
figure;plot(squeeze(errorspg(1,:,1)),squeeze(erroramp(:,:,4))./squeeze(erroramp(:,:,5)),'-','LineWidth',2);title('AMP')
xlabel('N');ylabel('relative error');%legend(num2str(matrixvectors(1)),num2str(matrixvectors(2)),num2str(matrixvectors(3)))%,num2str(iterations(4)),num2str(iterations(5)),num2str(iterations(6)),num2str(iterations(7)),num2str(iterations(8)),num2str(iterations(9)),num2str(iterations(10)))
axis tight

%%
figure('Name','Solution paths')
plot(infospg{end,end}.xNorm1,infospg{end,end}.rNorm2,infopqn{end,end}.xNorm1,infopqn{end,end}.rNorm2)%,infoamp{end,end}.xNorm1,infoamp{end,end}.rNorm2);hold on
scatter(infospg{end,end}.xNorm1,infospg{end,end}.rNorm2);
scatter(infopqn{end,end}.xNorm1,infopqn{end,end}.rNorm2);
scatter(infoamp{end,end}.xNorm1,infoamp{end,end}.rNorm2);hold off
xlabel('one-norm model')
ylabel('two-norm residual')
title('Solutions paths')
legend('SPGl1','PQNl1','AMP')
axis tight

%%
% figure(4);plot(squeeze(errorspg(1,:,2)),squeeze(errorspg(:,:,4))./squeeze(errorspg(:,:,5)),'-','LineWidth',2);title('SPGl1')
% xlabel('N');ylabel('relative error');legend(num2str(iterations(1)),num2str(iterations(2)),num2str(iterations(3)),num2str(iterations(4)),num2str(iterations(5)),num2str(iterations(6)),num2str(iterations(7)),num2str(iterations(8)),num2str(iterations(9)),num2str(iterations(10)))
% 
% figure(5);plot(squeeze(errorspg(1,:,1)),squeeze(erroramp(:,:,4))./squeeze(erroramp(:,:,5)),'-','LineWidth',2);title('AMP')
% xlabel('N');ylabel('relative error');legend(num2str(iterations(1)),num2str(iterations(2)),num2str(iterations(3)),num2str(iterations(4)),num2str(iterations(5)),num2str(iterations(6)),num2str(iterations(7)),num2str(iterations(8)),num2str(iterations(9)),num2str(iterations(10)))

%%
% figure;print -depsc2 Fig/varyingiterations_spg.eps
% figure;print -depsc2 Fig/varyingiterations_pqn.eps
% figure;print -depsc2 Fig/varyingiterations_amp.eps
% figure;print -depsc2 Fig/varyingiterations_path.eps