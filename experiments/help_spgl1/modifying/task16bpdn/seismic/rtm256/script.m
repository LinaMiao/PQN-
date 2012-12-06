%% add tools
addpath(genpath('/Volumes/Users/linamiao/Dropbox/extended_migration/data/Compass'))
addpath(genpath('/Volumes/Users/linamiao/Dropbox/extended_migration/functions'))
addpath(genpath('/users/slic/linamiao/Dropbox/extended_migration/data/Compass'))
addpath(genpath('/users/slic/linamiao/Dropbox/extended_migration/functions'))
cd ../../../../../../functions/
addpath(genpath(pwd))
cd ../experiments/help_spgl1/modifying/task16bpdn/seismic/
addpath(genpath('/users/slic/linamiao/Documents/Documents/Tools/Matlabtools/tristan/'))

% addpath(genpath('e:\research\Tools\tristan'))
% addpath(genpath('e:\research\Tools\spot-slim'))
% addpath(genpath('e:\research\Tools\pSPOT'))

%% generate data
    %% define model
    model.o = [0 0];
    model.d = [10 10];
    n1 = 32; n2 = 32;
    model.n = [n1 n2];
    model.nb = [30 30 0];
        
    t    = 0:.004:1; nt = length(t);
    freq = fkk(t); nf = length(freq);
    If   = [8:2:21];
    
    model.freq = freq(If);             nfreq = length(model.freq);
    model.zsrc = model.d(1);
    model.xsrc = linspace(0,(n2-1)*model.d(2),n2);  nsrc  = length(model.xsrc);
    model.zrec = model.d(1);
    model.xrec = linspace(0,(n2-1)*model.d(2),n2);  nrec  = length(model.xrec);
    model.f0   = 15;
    model.t0   = .01;
    
    xr = model.xrec; nr = length(xr);
    xs = model.xsrc; ns = length(xs);
    
    Q = speye(nsrc);
    
    z = 0:model.d(1):(n1-1)*model.d(1);nz = length(z);
    x = 0:model.d(2):(n2-1)*model.d(2);nx = length(x);
    [zz,xx] = ndgrid(z,x);
    
  
    %% define velocity
    v0 = 2000 + 0.*xx + 0.*zz;
    
    dv = 0.*xx;
   
    dv(zz >= 140) = 100;
    dv(zz >= 160) = 0;
    
    dv(xx <= 150) = 0;
    dv(xx >= 170) = 0;
    
%     dv(zz >= -100*cos(2*pi/(max(x)*1.5).*xx)+420) = 100;
%     dv(zz >= -100*cos(2*pi/(max(x)*1.5).*xx)+440) = 0;
    
   
%     dv(zz >= -100*cos(2*pi/(max(x)*1.5).*xx)+490) = 100;
%     dv(zz >= -100*cos(2*pi/(max(x)*1.5).*xx)+510) = 0;
    
    
%     dv(zz >= 100*sin(2*pi/(max(x)*1.7).*xx)+490) = 100;
%     dv(zz >= 100*sin(2*pi/(max(x)*1.7).*xx)+510) = 0;
%     
%     
%     dv(zz >= 100*sin(2*pi/(max(x)*1.5).*xx)+300) = 100;
%     dv(zz >= 100*sin(2*pi/(max(x)*1.5).*xx)+320) = 0;
%      
    v = v0 + dv;


    h = figure; imagesc(dv); colormap(gray);pbaspect([model.d(2),model.d(1),1]);title('dv');
    caxis([-1 1]*max(mean(dv)));
    saveas(h,'vel_pertb');
    
    m0 = 1e6./v0(:).^2;
    m = 1e6./v(:).^2;
    dm = m-m0;
    figure;imagesc(v);colormap(gray);pbaspect([model.d(2),model.d(1),1]);title('v')
    figure;imagesc(v0);colormap(gray);pbaspect([model.d(2),model.d(1),1]);title('v0')

    h = figure; imagesc(reshape(dm,model.n)); colormap(gray);pbaspect([model.d(2),model.d(1),1]);title('dm');
    caxis([-1 1]*abs(max(mean(reshape(dm,model.n)))));
    saveas(h,'m_pertb');
    %% data
    % generate data
    [~,J] = F(m0,Q,model);
    
    % linearize data
    Dobs = J*dm;
    Dobs = gather(Dobs);
  
    %save temp Dobs J
    %load temp
    % permute into f-x-x order
    cube = permute(reshape((Dobs),nrec,nsrc,length(model.freq)),[3 1 2]); 

    temp = zeros(nf,nrec,nsrc);
    temp(If,:,:) = cube([1:length(model.freq)],:,:);
    CUBE = temp;
    
%     % rtm
    rtm = J'*Dobs;
    rtm = reshape(rtm,nz,nr);
    h = figure;imagesc(real(rtm));colormap(gray);pbaspect([model.d(2),model.d(1),1]);title('rtm')
    saveas(h,'rtm')
% 
%     % lsqr
%     lsm = lsqr(J,Dobs,[],10);
%     lsm = reshape(lsm,nz,nr);
%     h = figure;imagesc(real(lsm));
%     saveas(h,'lsm')




%% migration

%C = opCurvelt(length(z),nr,6,16,1,'ME',0);
%C = opDFT2(nz,nr);
A = J;%*C';
b = Dobs;


% l1 migration
opts.iterations = 80;
opts.fid = fopen('log_freq_spg.txt', 'w'); 
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, b, 0, 1e-3, zeros(size(A,2),1), opts);
%m_spg = C'*x_spg;
m_spg = x_spg;
figure;imagesc(reshape(real(m_spg),nz,nr));
opts.fid = fopen('log_freq_pqn.txt', 'w'); 
sigma_ref = info_spg.rNorm;
corrections_list = [1 2 3 4 5];
for i = 3:length(corrections_list)
    corrections = corrections_list(i);
    [x_pqn(:,i),~,~,info_pqn(:,i)] = pqnl1_2(A, b, 0, 1e-3, zeros(size(A,2),1), opts,sigma_ref,corrections);
    %m_pqn(:,i) = C'*x_pqn(:,i);
    m_pqn(:,i) = x_pqn(:,i);
    nz = length(z);
    figure(i); subplot(2,1,1);imagesc(reshape(real(m_spg),nz,nr));
    subplot(2,1,2);imagesc(reshape(real(m_pqn(:,i)),nz,nr));
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


