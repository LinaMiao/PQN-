%% add tools
addpath(genpath('/Volumes/Users/linamiao/Dropbox/extended_migration/data/Compass'))
addpath(genpath('/Volumes/Users/linamiao/Dropbox/extended_migration/functions'))
addpath(genpath('/users/slic/linamiao/Dropbox/extended_migration/data/Compass'))
addpath(genpath('/users/slic/linamiao/Dropbox/extended_migration/functions'))
cd ../../../../../../functions/
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task16bpdn/seismic/rtm/
addpath(genpath('/users/slic/linamiao/Documents/Documents/Tools/Matlabtools/tristan/'))

% addpath(genpath('e:\research\Tools\tristan'))
% addpath(genpath('e:\research\Tools\spot-slim'))
% addpath(genpath('e:\research\Tools\pSPOT'))

%% generate data
    [v,o,d,n]=odnread('~/bg_vp.odn');
    v=reshape(v,n);
    v=v(1:2:end,1:2:end);
    v=v(1:101,1:101);

    v0=mean(v,2)*ones(1,101);
    v0=opSmooth(101,50)*v0;

    dv = v - v0;
    
    imagesc((v-v0)./v)
    m=1e6./v(:).^2;
    m0=1e6./v0(:).^2;
    dm=m-m0;

    imagesc(reshape(dm,101,101))

    
    [n1,n2] = size(v);

    %% define model
    model.o = [0 0];
    model.d = [10 10];
    model.n = [n1 n2];
    model.nb = [30 30 0];
        
    t    = 0:.004:1; nt = length(t);
    freq = fkk(t); nf = length(freq);
    If   = [8:1:26];
    
    model.freq = freq(If);             nfreq = length(model.freq);
    model.zsrc = model.d(1);
    model.xsrc = linspace(0,(n2-1)*model.d(2),n2);  nsrc  = length(model.xsrc);
    model.zrec = model.d(1);
    model.xrec = linspace(0,(n2-1)*model.d(2),n2);  nrec  = length(model.xrec);
    model.f0   = 10;
    model.t0   = .01;
    
    xr = model.xrec; nr = length(xr);
    xs = model.xsrc; ns = length(xs);
    
    Q = speye(nsrc);
    
    z = 0:model.d(1):(n1-1)*model.d(1);nz = length(z);
    x = 0:model.d(2):(n2-1)*model.d(2);nx = length(x);
    [zz,xx] = ndgrid(z,x);
    

    h = figure; imagesc(dm); colormap(gray);
    caxis([-1 1]*abs(max(max(reshape(dm,model.n)))));pbaspect([model.d(2),model.d(1),1]);
    saveas(h,'model_pertb');
    
  
    figure;imagesc(v);colormap(gray);pbaspect([model.d(2),model.d(1),1]);title('v')
    h = figure;imagesc(v0);colormap(gray);pbaspect([model.d(2),model.d(1),1]);title('v0')
    saveas(h,'v0');

    
    
    %% data
    % generate data
    [~,J] = F(m0,Q,model);
    
    % linearize data
    Dobs = J*dm;
    Dobs = gather(Dobs);
  
    %save temp Dobs J
    %load temp
    
    % rtm
    rtm = J'*Dobs;
    rtm = reshape(rtm,nz,nr);
    h = figure;imagesc(real(rtm));colormap(gray);pbaspect([model.d(2),model.d(1),1]);title('rtm')
    caxis([-1 1]*abs(max(max(reshape(dm,model.n)))));pbaspect([model.d(2),model.d(1),1]);
    saveas(h,'rtm')

    % lsqr
    lsm = lsqr(J,Dobs,[],10);
    lsm = reshape(lsm,nz,nr);
    h = figure;imagesc(real(lsm));
    saveas(h,'lsm')




%% migration

C = opCurvelet(length(z),nr,6,16,1,'ME',0);

A = J*C';
b = Dobs;


% l1 migration
opts.iterations = 100;
opts.fid = fopen('log_freq_spg.txt', 'w'); 
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, b, 0, 1e-3, zeros(size(A,2),1), opts);
m_spg = C'*x_spg;
h = figure;imagesc(reshape(real(m_spg),nz,nr));    caxis([-1 1]*abs(max(max(reshape(dm,model.n)))));pbaspect([model.d(2),model.d(1),1]);
saveas(h,'mspg')
opts.fid = fopen('log_freq_pqn.txt', 'w'); 
sigma_ref = info_spg.rNorm;
corrections_list = [5 10 15 20 25 30];
for i = 1:length(corrections_list)
    corrections = corrections_list(i);
    [x_pqn(:,i),~,~,info_pqn(:,i)] = pqnl1_2(A, b, 0, 1e-3, zeros(size(A,2),1), opts,sigma_ref,corrections);
    m_pqn(:,i) = C'*x_pqn(:,i);
    nz = length(z);
    h = figure(i); subplot(2,1,1);imagesc(reshape(real(m_spg),nz,nr));
    caxis([-1 1]*abs(max(max(reshape(dm,model.n)))));pbaspect([model.d(2),model.d(1),1]);
    subplot(2,1,2);imagesc(reshape(real(m_pqn(:,i)),nz,nr));
    caxis([-1 1]*abs(max(max(reshape(dm,model.n)))));pbaspect([model.d(2),model.d(1),1]);
    saveas(h,'spg vs pqn')
    save info info_spg info_pqn m_spg m_pqn
end
    

% show solution path for the last one
h = figure('Name','Solution paths');
plot(info_spg.xNorm1,info_spg.rNorm2,info_pqn(:,i).xNorm1,info_pqn(:,i).rNorm2);hold on
scatter(info_spg.xNorm1,info_spg.rNorm2);
scatter(info_pqn(:,i).xNorm1,info_pqn(:,i).rNorm2);hold off
legend('SPGL1','PQNl1')
axis tight
saveas(h,'solution path freq.jpg')


