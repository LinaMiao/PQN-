%% add tools
cd ./functions
addpath(genpath(pwd))
cd ../../../../../../../functions/
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task17mmv/seismic/l1migration/
addpath(genpath('/users/slic/linamiao/Documents/Documents/Tools/Matlabtools/tristan/'))

% addpath(genpath('e:\research\Tools\tristan'))
% addpath(genpath('e:\research\Tools\spot-slim'))
% addpath(genpath('e:\research\Tools\pSPOT'))

%% generate data
 %% define model
    model.o = [0 0];
    model.d = [10 10];
    model.n = [64 64];
    model.nb = [10 10 0];
    model.t0 = -.08;
    
    t    = 0:.004:1; nt = length(t);
    freq = fkk(t); nf = length(freq);
    If   = [5:1:15];
    
    model.freq = freq(If);             nfreq = length(model.freq);
    model.zsrc = 10;
    model.xsrc = linspace(0,630,16);  nsrc  = length(model.xsrc);
    model.zrec = 10;
    model.xrec = linspace(0,630,16);  nrec  = length(model.xrec);
    model.f0   = 10;
    model.t0   = -.08;
    
    xr = model.xrec; nr = length(xr);
    xs = model.xsrc; ns = length(xs);
    
    Q = speye(nsrc);
    
    z = 0:10:630; nz = length(z);
    x = 0:10:630; nx = length(x);
    [zz,xx] = ndgrid(z,x);
    
    % background velocity [m/s]
    v0 = 2000 + 0.*zz;
    dv = 0*xx;
    dv(zz >= 300) = 100;
    dv(zz >= 360) = 0;
    figure; imagesc(dv);
    v = v0+dv;
    m0 = 1e6./v0(:).^2;
    m = 1e6./v(:).^2;
    dm = m-m0;
    

    
    

    %% data
    % generate data
    [~,J] = F(m0,Q,model);
    
    % linearize data
    Dobs = J*dm;
    Dobs = gather(Dobs);
  
    save temp Dobs
    load temp
    % permute into f-x-x order
    cube = permute(reshape((Dobs),nrec,nsrc,length(model.freq)),[3 1 2]); 

    temp = zeros(nf,nrec,nsrc);
    temp(If,:,:) = cube([1:length(model.freq)],:,:);
    CUBE = temp;

    % window data
    %CUBE = CUBE(:,1:5:end,:);
    xr = xs;
    nr = length(xr);
    

    % what the data look like in time domain
    Dt = ifft(CUBE,[],1)*sqrt(size(CUBE,1));
    
    figure;imagesc(squeeze(real(Dt(:,:,8))));




%% migration

% %% least square migration
temp = reshape(dm,model.n);
temp = temp(:,1:4:end);
temp = vec(temp);
% 

% in time domain
dim = size(Dt);
At = opDSR_mig_time(dim,t,xr,xs,z,2000*ones(length(z),1),0);
Bt = At'*temp;
[Xt] = lsqr(At',Bt,[],1); 
lsxt = reshape(Xt,model.n(1),length(xr));
figure;imagesc(real(lsxt));title('time domain least square migration')
% adjoint result
adxt = At * Bt;
adxt = reshape(adxt,model.n(1),length(xr));
figure;imagesc(real(adxt));title('adjoint of operator')


%% simutaneous experiment with randomized source superposition
% source superposition
ncs = 5;
M   = opGaussian(ns,ncs);
RM  = opKron(opDirac(nt),opDirac(nr),M');

simD = RM*vec(Dt);
simD = reshape(simD,nt,nr,ncs);



%% extended image gather
dim = [nt,nr,ns];
option.freq = 0;
K = opDSR(dim,t,xr,xs,z,v,option);
M = K*vec(Dt); 

ny = nr; % number of midpoints
np = numel(M)/(ny*nz); % number of wavenumbers
M = reshape(M,ny,nz,np); 
figure;
for i = 1:ny
    imageHandle = subplot(1,ny,i); imagesc(real(squeeze(M(i,:,:))));colormap(gray);axis off;
    pHandle = get(imageHandle, 'pos');
    pHandle(3) = pHandle(3) + 0.005;
    set(imageHandle, 'pos', pHandle);
end


%% joint sparsity
% sparse transform
ip = 600;
mm = reshape(M(:,:,ip),ny,nz);
norm(vec(mm),1)/length(vec(mm))

C1 = opCurvelet(nz,ny,6,16,1,'ME',0);
aa1 = C1 * vec(mm);
norm(aa1,1)/length(aa1)

C2 = opFFT2(size(mm));
aa2 = C2 * vec(mm);
norm(aa2,1)/length(aa2)

h = figure;
subplot(3,1,1); plot(sort(abs(vec(mm)),'descend')); title('original');
subplot(3,1,2); plot(sort(abs(vec(aa1)),'descend')); title('curvelet')
subplot(3,1,3); plot(sort(abs(vec(aa2)),'descend')); title('2DFFT')
saveas(h,'sparsify transform.jpg')

%%
C = C2;
S = opKron(opDirac(np),C1);
A = RM*K'*S';
B = vec(simD);

groups = np;

opts.fid = fopen('log_log.txt', 'w'); 
%opts.iterations = 1;
opts.project     = @(x,weight,tau) NormL12_project(groups,x,weight,tau);
opts.primal_norm = @(x,weight    ) NormL12_primal(groups,x,weight);
opts.dual_norm   = @(x,weight    ) NormL12_dual(groups,x,weight);


sigma = 1e-3;
opts.iterations = 1;
opts.fid = fopen('spg','w');
[x_spg,r_spg,g_spg,info_spg] = spgl1(A,B(:),0,sigma,zeros(size(A,2),1),opts);
opts.fid = fopen('pqn','w');
%sigma_ref = info_spg.rNorm;
sigma_ref = 1e2;
opts.decTol = 1e-4;
opts.iterations = 20;
[x_pqn1,r_pqn1,g_pqn1,info_pqn1] = pqnl1_2(A,B(:),0,sigma,zeros(size(A,2),1),opts,sigma_ref);
info_spg
info_pqn

% show optimiztion result
h = figure('Name','Solution paths');
plot(info_spg.xNorm1,info_spg.rNorm2,info_pqn.xNorm1,info_pqn.rNorm2);hold on
scatter(info_spg.xNorm1,info_spg.rNorm2);
scatter(info_pqn1.xNorm1,info_pqn1.rNorm2);hold off
legend('SPGL1','PQNl1')
axis tight
m_spg = S'*x_spg;
m_pqn = S'*x_pqn1;
save info_freq info_spg info_pqn m_spg m_pqn
saveas(h,'solution path freq.jpg')


%% Coommon p image gather->Image
op = opImgather([nz,nr,ns],z,xr,xs,t);
exIm_spg = op' * m_spg;
exIm_spg = reshape(exIm_spg,nz,nr,ns);
Im_spg = zeros(nz,nr);
for i = 1:nz
    Im_spg(i,:) = diag(squeeze(exIm_spg(i,:,:)));
end

exIm_pqn = op' * m_pqn;
exIm_pqn = reshape(exIm_pqn,nz,nr,ns);
Im_pqn = zeros(nz,nr);
for i = 1:nz
    Im_pqn(i,:) = diag(squeeze(exIm_pqn(i,:,:)));
end


h = figure;
subplot(1,2,1); imagesc(real(Im_spg)); colormap(gray);title('spg')
subplot(1,2,2); imagesc(real(Im_pqn)); colormap(gray);title('pqn')
saveas(h,'joint sparsity image.jpg')


