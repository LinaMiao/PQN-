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
    load compass_velocity 
    compass = Data;
    clear Data;
    v = compass(end-180:1:end-80,400:1:500);
    temp = smooth(v);
    temp = reshape(temp,size(v));
    avg = mean(temp,2);
    v0 = zeros(size(v));
    for i = 1:size(v,2)
        v0(:,i) = avg;
    end
    clear avg
    dv = v-v0;
    
    [n1,n2] = size(v);

    %% define model
    model.o = [0 0];
    model.d = [6 25];
    model.n = [n1 n2];
    model.nb = [30 30 0];
        
    t    = 0:.004:1; nt = length(t);
    freq = fkk(t); nf = length(freq);
    If   = [8:1:20];
    
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
    

    h = figure; imagesc(dv); colormap(gray);pbaspect([model.d(2),model.d(1),1]);title('dv');
    saveas(h,'vel_pertb');
    
    m0 = 1e6./v0(:).^2;
    m = 1e6./v(:).^2;
    dm = m-m0;
    figure;imagesc(v);colormap(gray);pbaspect([model.d(2),model.d(1),1]);title('v')
    h = figure;imagesc(v0);colormap(gray);pbaspect([model.d(2),model.d(1),1]);title('v0')
    saveas(h,'v0');

    
    

    %% data
    % generate data
%     [~,J] = F(m0,Q,model);
%     
%     % linearize data
%     Dobs = J*dm;
%     Dobs = gather(Dobs);
%   
%     save temp Dobs J
    load temp
    % permute into f-x-x order
    cube = permute(reshape((Dobs),nrec,nsrc,length(model.freq)),[3 1 2]); 

    temp = zeros(nf,nrec,nsrc);
    temp(If,:,:) = cube([1:length(model.freq)],:,:);
    CUBE = temp;






%% migration

C = opCurvelet(length(z),nr,6,16,1,'ME',0);

A = J*C';
b = Dobs;


% l1 migration
opts.iterations = 30;
opts.fid = fopen('log_freq_spg.txt', 'w'); 
[x_spg,r_spg,g_spg,info_spg] = spgl1(A, b, 0, 1e-3, zeros(size(A,2),1), opts);
m_spg = C'*x_spg;
opts.fid = fopen('log_freq_pqn.txt', 'w'); 
sigma_ref = info_spg.rNorm;
corrections_list = [5 10 15 20 25 30];
for i = 1:length(corrections_list)
    corrections = corrections_list(i);
    [x_pqn(:,i),~,~,info_pqn(:,i)] = pqnl1_2(A, b, 0, 1e-3, zeros(size(A,2),1), opts,sigma_ref,corrections);
    m_pqn(:,i) = C'*x_pqn(:,i);
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


