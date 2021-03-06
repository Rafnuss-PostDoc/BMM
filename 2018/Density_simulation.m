%% Simulation Map
% Script to generate simulated map of bird migration

%% Clear and load data
clear all;
% load('./data/dc_corr.mat');
load('./data/Density_inference.mat');
load('data/Density_estimationMap','g')
load coastlines;
addpath('./functions/');


%% FOr simulating more realization for a few night only.
% Amplitude calculation is the same (i.e, we simulate all day), but
% residuals are only computed for the day needed
% Uncomment all text in the code with `nightly_cmt`
% nightly_i_t=52:55; %nightly_cmt
% sim_nightly = ismember(g.day_id,nightly_i_t); %nightly_cmt

%% 1. multitude

% 1.1. Creation of the grid and path
[LAT,LON,AT]= ndgrid(1:g.nlat,1:g.nlon,1:g.nat);
Path=Inf*ones(g.nlat,g.nlon,g.nat);
Path(repmat(g.latlonmask,1,1,g.nat))=nan;

slat = 1:ceil(log(g.nlat+1)/log(2));
slon = 1:ceil(log(g.nlon+1)/log(2));
sat = 1:ceil(log(g.nat+1)/log(2));
sn = max([numel(slat), numel(slon) numel(sat)]);
nb = nan(sn,1);
start = zeros(sn+1,1);
path = nan(sum(isnan(Path(:))),1);
ds = 2.^(sn-1:-1:0);
for i_scale = 1:sn
    [LAT_s,LON_s,AT_s] = ndgrid(1:ds(i_scale):g.nlat,1:ds(i_scale):g.nlon,1:ds(i_scale):g.nat); % matrix coordinate
    id = find(isnan(Path(:)));
    id = id(ismember([LAT(id) LON(id) AT(id)], [LAT_s(:) LON_s(:) AT_s(:)], 'rows'));
    nb(i_scale) = numel(id);
    start(i_scale+1) = start(i_scale)+nb(i_scale);
    path( start(i_scale)+(1:nb(i_scale)) ) = id(randperm(nb(i_scale)));
    Path(path( start(i_scale)+(1:nb(i_scale)) )) = start(i_scale)+(1:nb(i_scale));
end

Pathll=reshape(Path(repmat(g.latlonmask,1,1,g.nat)),g.nlm,g.nat);
[~,pathll] = sort(Pathll(:));

% 1.2 Ditance Matrix
Ddist_gg = squareform(pdist([g.lat2D(g.latlonmask) g.lon2D(g.latlonmask)], @lldistkm));
Ddist_gr = pdist2([data.lat(multi.R(:)) data.lon(multi.R(:))], [g.lat2D(g.latlonmask) g.lon2D(g.latlonmask)], @lldistkm);
Ddist_gr = Ddist_gr(~multi.isnan(:),:);
Dtime_gg = squareform(pdist(datenum(g.day)));
Dtime_gr = pdist2(multi.T(:), datenum(g.day));
Dtime_gr = Dtime_gr(~multi.isnan,:);
Dtime_rr = squareform(pdist(multi.T(:)));
Dtime_rr = Dtime_rr(~multi.isnan(:),~multi.isnan(:));
Ddist_rr = squareform(pdist([data.lat(multi.R(:)) data.lon(multi.R(:))], @lldistkm));
Ddist_rr = Ddist_rr(~multi.isnan(:),~multi.isnan(:));

C_rr(:,:,1) = multi.cov.Cn(Ddist_rr,Dtime_rr,multi.cov.parm(:,1)) + diag(multi.M_var(~multi.isnan));
C_rr(:,:,2) = multi.cov.Cn(Ddist_rr,Dtime_rr,multi.cov.parm(:,2)) + diag(multi.M_var(~multi.isnan));
C_rr(:,:,3) = multi.cov.Cn(Ddist_rr,Dtime_rr,multi.cov.parm(:,3)) + diag(multi.M_var(~multi.isnan));
C_rr(:,:,4) = multi.cov.Cn(Ddist_rr,Dtime_rr,multi.cov.parm(:,4)) + diag(multi.M_var(~multi.isnan));

% 1.3. Initizialization of the kriging weights and variance error,Computation of neighbours, lambda weight and S
neighg_nb=100;
neighr_nb=5;
NEIGHG = nan(neighg_nb,numel(path),'single');
NEIGHR = nan(neighr_nb,numel(path),'single');
LAMBDA = nan(neighg_nb+neighr_nb,numel(path),'single');
S = nan(numel(path),1,'single');
tmp = repmat((1:g.nlm)',1,g.nat);
LL_i=tmp(pathll);
tmp = repmat(1:numel(g.day),g.nlm,1);
TT_i=tmp(pathll);

wradius=3;

for i_pt = 1:numel(path)
    
    % Get season id
    i_b = g.day_b(TT_i(i_pt));
    
    % 1. Neighbor from grid
    % Find in the index of the path such that the spatio-temporal distance
    % is within wradius. Then, only select the already simulated value
    % (path value less than currently simulated)
    neighg=find(Pathll<i_pt & bsxfun(@and,Ddist_gg(:,LL_i(i_pt))<multi.cov.parm(4,i_b)*wradius  , Dtime_gg(:,TT_i(i_pt))'<multi.cov.parm(5,i_b)*wradius));
    Cgp = multi.cov.C(Ddist_gg(LL_i(Pathll(neighg)),LL_i(i_pt)), Dtime_gg(TT_i(Pathll(neighg)),TT_i(i_pt)),multi.cov.parm(:,i_b));
    [Cgp,tmp]=maxk(Cgp,neighg_nb);
    neighg=neighg(tmp);
    Cgg = multi.cov.C(Ddist_gg(LL_i(Pathll(neighg)),LL_i(Pathll(neighg))), Dtime_gg(TT_i(Pathll(neighg)),TT_i(Pathll(neighg))), multi.cov.parm(:,i_b));
    
    % 2. Find the radar neigh
    neighr=find(Ddist_gr(:,LL_i(i_pt))<multi.cov.parm(4,i_b)*wradius & Dtime_gr(:,TT_i(i_pt))<multi.cov.parm(5,i_b)*wradius);
    Crp = multi.cov.C(Ddist_gr(neighr,LL_i(i_pt)), Dtime_gr(neighr,TT_i(i_pt)),multi.cov.parm(:,i_b));
    [Crp,tmp] = maxk(Crp,neighr_nb);
    neighr=neighr(tmp);
    Crr = C_rr(neighr,neighr,i_b);
    
    Crg = multi.cov.C(Ddist_gr(neighr,LL_i(Pathll(neighg))), Dtime_gr(neighr,TT_i(Pathll(neighg))),multi.cov.parm(:,i_b));
    
    l =  [Cgg Crg'; Crg Crr] \ [Cgp; Crp];
    
    NEIGHG(:,i_pt) = [neighg ; nan(neighg_nb-numel(neighg),1)];
    NEIGHR(:,i_pt) = [neighr ; nan(neighr_nb-numel(neighr),1)];
    LAMBDA(:,i_pt) = [l ; nan(neighr_nb+neighg_nb-numel(l),1)]';
    
    S(i_pt) = sqrt( multi.cov.C(0,0, multi.cov.parm(:,i_b)) - l'*[Cgp; Crp]);
    assert(isreal(S(i_pt)))
end


% save('D:\Guests\rafnuss_at_gmail_com\tmp_BMM2018\Density_simulationMap_multi','NEIGHG','NEIGHR','LAMBDA','pathll','S')
% load('D:\Guests\rafnuss_at_gmail_com\tmp_BMM2018\Density_simulationMap_multi')

%% 2. Amplitude: Simulation of realization


nb_real = 100;
real_M_mn_ll = nan(g.nlm,g.nat,nb_real);
multi_M_mn = multi.M_mn(~multi.isnan);

for i_real=1:nb_real
    real_multi_tmp=nan(g.nlm,g.nat);
    rng('shuffle');
    U=randn(g.nlm,g.nat);
    
    for i_pt = 1:numel(pathll)
        ng = ~isnan(NEIGHG(:,i_pt));
        nr = ~isnan(NEIGHR(:,i_pt));
        nl = ~isnan(LAMBDA(:,i_pt));
        real_multi_tmp(pathll(i_pt)) = LAMBDA(nl,i_pt)'*[real_multi_tmp(NEIGHG(ng,i_pt)) ; multi_M_mn(NEIGHR(nr,i_pt))] + U(i_pt)*S(i_pt);
        assert(~isnan(real_multi_tmp(pathll(i_pt))))
        assert(isreal(real_multi_tmp(pathll(i_pt))))
    end
    real_M_mn_ll(:,:,i_real) = real_multi_tmp;
end

real_An = nan(g.nlat,g.nlon,g.nat,nb_real);
real_An(repmat(g.latlonmask,1,1,g.nat,nb_real)) = real_M_mn_ll;
real_A = nan(g.nlat,g.nlon,g.nat,nb_real);

for i_b=1:numel(data.block.date)
    id = g.day_b==i_b;
    real_A(:,:,id,:) = real_An(:,:,id,:) + multi.f_trend(permute((repmat(datenum(g.day(id)),1,g.nlat,g.nlon,nb_real)),[2 3 1 4]), ...
        repmat(g.lat2D,1,1,sum(id),nb_real), repmat(g.lon2D,1,1,sum(id),nb_real), multi.beta(:,i_b));
end

%save('D:\Guests\rafnuss_at_gmail_com\tmp_BMM2018\real_A','real_A' ,'-v7.3')

clear real_An real_M_mn_ll multi_M_mn U

figure('position',[0 0 1000 600]);
for i_t = 1:g.nat-1
    subplot(4,6,i_t); hold on; worldmap([min(data.lat) max(data.lat)], [min(data.lon) max(data.lon)]);
    surfm(g.lat2D,g.lon2D,real_A(:,:,i_t))
    scatterm(data.lat,data.lon,[],multi.M_m(i_t,:),'filled','MarkerEdgeColor','k')
    plotm(coastlat, coastlon,'k')
end

% fig=figure(2);
% filename='data/Density_simulationMap_amplitude';
% Frame(g.nat-1) = struct('cdata',[],'colormap',[]);
% worldmap([min(data.lat) max(data.lat)], [min(data.lon) max(data.lon)]);
% hcoast=plotm(coastlat, coastlon,'k'); caxis([min(multi.A(:)) max(multi.A(:))]);
% for i_t = 1:g.nat-1
%
%     hsurf=surfm(g.lat2D,g.lon2D,real_A(:,:,i_t));
%     hscat=scatterm(data.lat,data.lon,[],multi.A(:,i_t),'filled','MarkerEdgeColor','k');
%     uistack(hcoast,'top')
%
%     title(datestr(g.atime(i_t))); drawnow;
%
%     Frame(i_t) = getframe(fig);
%     [imind,cm] = rgb2ind(frame2im(Frame(i_t)),256);
%     if i_t == 1
%         imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append');
%     end
%     delete(hsurf); delete(hscat);
% end

















%% 2. Residu:

% 2.1. Residu: Generation of the path
mask_total = repmat(g.latlonmask,1,1,g.nt);
mask_total(mask_total) = g.mask_day(:);

Pathll=cell(1,g.nat);
pathll=cell(1,g.nat);
ndt=cell(1,g.nat);

slat = 1:ceil(log(g.nlat+1)/log(2));
slon = 1:ceil(log(g.nlon+1)/log(2));

tic  %140sec not-parr because of ismember (lin 220)
for i_t=1:g.nat
    ndt{i_t}=sum(g.day_id==i_t);
    [LAT,LON,DT]= ndgrid(1:g.nlat,1:g.nlon,1:ndt{i_t});
    
    Path=Inf*ones(g.nlat,g.nlon,ndt{i_t});
    Path(mask_total(:,:,g.day_id==i_t))=nan;
    path = nan(sum(isnan(Path(:))),1);
    
    sat = 1:ceil(log(ndt{i_t}+1)/log(2));
    sn = max([numel(slat), numel(slon) numel(sat)]);
    nb = nan(sn,1);
    start = zeros(sn+1,1);
    ds = 2.^(sn-1:-1:0);
    for i_scale = 1:sn
        [LAT_s,LON_s,DT_s] = ndgrid(1:ds(i_scale):g.nlat,1:ds(i_scale):g.nlon,1:ds(i_scale):ndt{i_t}); % matrix coordinate
        id = find(isnan(Path(:)));
        id = id(ismember([LAT(id) LON(id) DT(id)], [LAT_s(:) LON_s(:) DT_s(:)], 'rows'));
        nb(i_scale) = numel(id);
        start(i_scale+1) = start(i_scale)+nb(i_scale);
        path( start(i_scale)+(1:nb(i_scale)) ) = id(randperm(nb(i_scale)));
        Path(path( start(i_scale)+(1:nb(i_scale)) )) = start(i_scale)+(1:nb(i_scale));
    end
    
    Pathll{i_t}=reshape(Path(repmat(g.latlonmask,1,1,ndt{i_t})), g.nlm, ndt{i_t});
    [~,pathll{i_t}] = sort(Pathll{i_t}(:));
    pathll{i_t}=pathll{i_t}(~isinf(Pathll{i_t}(pathll{i_t})));
    
    Pathll{i_t}= int64(Pathll{i_t});
    pathll{i_t}= int64(pathll{i_t});
end
toc

% Distance Matrix common to all day
Ddist_gg = squareform(pdist([g.lat2D(g.latlonmask) g.lon2D(g.latlonmask)], @lldistkm));
Ddist_gr = pdist2([data.lat data.lon], [g.lat2D(g.latlonmask) g.lon2D(g.latlonmask)], @lldistkm);

% Reconstruct some intra variable not saved in mat file to save space.
intra_isnan = isnan(intra.I_m);
tmp=repmat(1:data.nrad,data.ntime,1);
intra_radar = tmp(~intra_isnan);
tmp=repmat(data.day(data.day_id),1,data.nrad);
intra_day = tmp(~intra_isnan);
tmp=repmat(datenum(data.time),1,data.nrad);
intra_time = tmp(~intra_isnan);

%%
% 2.2 Residu: Computaiton of the weight neigh and S

NEIGHG=cell(g.nat,1);
NEIGHR=cell(g.nat,1);
LAMBDA=cell(g.nat,1);
S=cell(g.nat,1);
neighday=cell(g.nat,1);
intra_cov_C=cell(g.nat,1);
sim=cell(g.nat,1);

wradius=4;
neighg_nb=100;
neighr_nb=50;


% for i_t=nightly_i_t %nightly_cmt
for i_t=1:g.nat
    i_b = g.day_b(i_t);
    intra_cov_parm = intra.cov.parm(:,i_b);
    
    % Sub-selection of the point to simulate during this day
    sim{i_t} = find(g.day_id==i_t);
    
    % Sub-selection of the radar data during this day
    neighday{i_t}  = find(intra_day==datenum(g.day(i_t)));
    Dtime = squareform(pdist(intra_time(neighday{i_t})));
    Ddist = squareform(pdist([data.lat(intra_radar(neighday{i_t})) data.lon(intra_radar(neighday{i_t}))], @lldistkm));
    intra_cov_C_i_d = intra.cov.Cn(Ddist,Dtime,intra_cov_parm);
    
    % Distance matrix for the day
    Dtime_gg = squareform(pdist(datenum(g.time(sim{i_t}))));
    Dtime_gr = pdist2(intra_time(neighday{i_t}), datenum(g.time(sim{i_t})));
    
    % 4. Initizialization of the kriging weights and variance error
    NEIGHG_tmp = nan(neighg_nb,numel(pathll{i_t}));
    NEIGHR_tmp = nan(neighr_nb,numel(pathll{i_t}));
    LAMBDA_tmp = nan(neighg_nb+neighr_nb,numel(pathll{i_t}));
    S_tmp = nan(numel(pathll{i_t}),1);
    tmp = repmat((1:g.nlm)',1,ndt{i_t});
    LL_i=tmp(pathll{i_t});
    tmp = repmat(1:ndt{i_t},g.nlm,1);
    TT_i=tmp(pathll{i_t});
    
    Pathll_i_d = Pathll{i_t};
    
    % 5 Loop of scale for multi-grid path
    tic
    parfor i_pt = 1:numel(pathll{i_t})
        
        % 1. Neighbor from grid
        % Find in the index of the path such that the spatio-temporal distance
        % is within wradius. Then, only select the already simulated value
        % (path value less than currently simulated)
        neighg = find( Pathll_i_d<i_pt & bsxfun(@and,Ddist_gg(:,LL_i(i_pt))<intra_cov_parm(4)*wradius  , Dtime_gg(:,TT_i(i_pt))'<intra_cov_parm(5)*wradius));
        Cgp = intra.cov.C(Ddist_gg(LL_i(Pathll_i_d(neighg)),LL_i(i_pt)), Dtime_gg(TT_i(Pathll_i_d(neighg)), TT_i(i_pt)), intra_cov_parm);
        [Cgp,tmp]=maxk(Cgp,neighg_nb);
        neighg=neighg(tmp);
        
        Cgg = intra.cov.C(Ddist_gg(LL_i(Pathll_i_d(neighg)),LL_i(Pathll_i_d(neighg))), Dtime_gg(TT_i(Pathll_i_d(neighg)),TT_i(Pathll_i_d(neighg))), intra_cov_parm);
        
        % 2. Find the radar neigh
        neighr=find( Ddist_gr(intra_radar(neighday{i_t}),LL_i(i_pt))<intra_cov_parm(4)*wradius & Dtime_gr(:,TT_i(i_pt))<intra_cov_parm(5)*wradius);
        Crp = intra.cov.C(Ddist_gr(intra_radar(neighday{i_t}(neighr)),LL_i(i_pt)), Dtime_gr(neighr,TT_i(i_pt)), intra_cov_parm);
        [Crp,tmp]=maxk(Crp,neighr_nb);
        neighr=neighr(tmp);
        Crr = intra_cov_C_i_d(neighr,neighr);
        
        % 3. cross covariance
        Crg = intra.cov.C(Ddist_gr(intra_radar(neighday{i_t}(neighr)),LL_i(Pathll_i_d(neighg))), Dtime_gr(neighr,TT_i(Pathll_i_d(neighg))), intra_cov_parm);
        
        l =  [Cgg Crg'; Crg Crr] \  [Cgp; Crp];
        
        NEIGHG_tmp(:,i_pt) = [neighg ; nan(neighg_nb-numel(neighg),1)];
        NEIGHR_tmp(:,i_pt) = [neighr ; nan(neighr_nb-numel(neighr),1)];
        LAMBDA_tmp(:,i_pt) = [l ; nan(neighr_nb+neighg_nb-numel(l),1)]';
        
        S_tmp(i_pt) = sqrt( intra.cov.C(0,0,intra_cov_parm) - l'*[Cgp; Crp]);
        
    end
    
    NEIGHG{i_t} = single(NEIGHG_tmp);
    NEIGHR{i_t} = single(NEIGHR_tmp);
    LAMBDA{i_t} = single(LAMBDA_tmp);
    S{i_t} = single(S_tmp);
    
    %     dlmwrite(['data/sim_intra_NEIGHG_' num2str(i_d)],NEIGHG{i_d})
    %     dlmwrite(['data/sim_intra_NEIGHR_' num2str(i_d)],NEIGHR{i_d})
    %     dlmwrite(['data/sim_intra_LAMBDA_' num2str(i_d)],LAMBDA{i_d})
    
    toc
    
end


% save('D:\Guests\rafnuss_at_gmail_com\tmp_BMM2018\Density_simulationMap_residu','NEIGHG','NEIGHR','LAMBDA','pathll','S','neighday','sim','-v7.3')
% load('D:\Guests\rafnuss_at_gmail_com\tmp_BMM2018\Density_simulationMap_residu')


%% 2.2 Residu: Generation of realizations
%load('data/Density_simulationMap_residu')

[latmask, lonmask]=ind2sub([g.nlat g.nlon],find(g.latlonmask));

intra_I_mn = intra.I_mn;
nb_real = 100;
real_In = nan(g.nlat,g.nlon,g.nt,nb_real,'single');
% real_In = nan(g.nlat,g.nlon,sum(sim_nightly),nb_real,'single'); %nightly_cmt

parpool(40)

for i_real=1:nb_real
    tic
    i_real
    real_intra_ll=cell(g.nat,1);
    rng('shuffle');
    % U=randn(g.nlm,g.nt);
    g_nlm = g.nlm;
    
    parfor i_d=1:g.nat
        
        tmp = nan(g_nlm,numel(sim{i_d}));
        for i_pt = 1:numel(pathll{i_d})
            ng = ~isnan(NEIGHG{i_d}(:,i_pt));
            nr = ~isnan(NEIGHR{i_d}(:,i_pt));
            nl = ~isnan(LAMBDA{i_d}(:,i_pt));
            tmp(pathll{i_d}(i_pt)) = LAMBDA{i_d}(nl,i_pt)'*[tmp(NEIGHG{i_d}(ng,i_pt)) ; intra_I_mn(neighday{i_d}(NEIGHR{i_d}(nr,i_pt)))] + randn*S{i_d}(i_pt);
            %assert(~isnan(tmp(pathll{i_d}(i_pt))))
            %assert(isreal(tmp(pathll{i_d}(i_pt))))
        end
        real_intra_ll{i_d} = tmp;
    end
    
    tmp2 = nan(g.nlm,g.nt);
    for i_d=1:g.nat
        tmp2(:,sim{i_d}) = real_intra_ll{i_d};
    end
    
    tmp3 = nan(g.nlat,g.nlon,g.nt);
    tmp3(repmat(g.latlonmask,1,1,g.nt))= tmp2;
    
    
    
    % tmp2 = nan(g.nlat,g.nlon,sum(sim_nightly)); %nightly_cmt
    % tmp2(repmat(g.latlonmask,1,1,sum(sim_nightly)))= real_intra_ll(:,sim_nightly); %nightly_cmt
    
    real_In(:,:,:,i_real) = tmp3;
    toc
end


real_I = nan(g.nlat,g.nlon,g.nt,nb_real,'single');
% real_I = nan(g.nlat,g.nlon,sum(sim_nightly),nb_real,'single');%nightly_cmt

g_NNT = nan(g.nlat,g.nlon,g.nt);
g_NNT(repmat(g.latlonmask,1,1,g.nt)) = g.NNT;

% i_b=1;%nightly_cmt

for i_b=1:numel(data.block.date)
    
    id = g.time_b==i_b;
    % id = sim_nightly;%nightly_cmt
    
    tmp = sqrt(intra.f_trend_v(g_NNT(:,:,id),intra.beta_v(:,i_b)));
    tmp2 = intra.f_trend(g_NNT(:,:,id),intra.beta(:,i_b));
    
    for i_real=1:nb_real
        real_I(:,:,id,i_real) = real_In(:,:,id,i_real) .* tmp + tmp2;
    end
    %real_I(:,:,id,:) = real_In(:,:,id,:) .* repmat(tmp,1,1,1,nb_real) + repmat(tmp2,1,1,1,nb_real);
    % real_I = real_In .* repmat(tmp,1,1,1,nb_real) + repmat(tmp2,1,1,1,nb_real); %nightly_cmt
end



%%


i_real=1;
h=figure(2);
worldmap([min(data.lat) max(data.lat)], [min(data.lon) max(data.lon)]);
filename='data/Density_simulationMap_residu';
mask_fullday = find(~reshape(all(all(isnan(real_I(:,:,:,i_real)),1),2),g.nt,1));
Frame(numel(mask_fullday)-1) = struct('cdata',[],'colormap',[]);
hcoast=plotm(coastlat, coastlon,'k'); caxis([-3 3]); hscat=[];
for i_t = 2:numel(mask_fullday)
    
    hsurf=surfm(g.lat2D,g.lon2D,real_In(:,:,mask_fullday(i_t),i_real));
    
    %     gtit = datenum(g.time(mask_fullday(i_t)+([-1 0 1])));
    %     id = find(mean(gtit(1:2))<data.time & mean(gtit(2:3))>data.time);
    %     if sum(id)>0
    %         G = findgroups(data.dateradar(id));
    %         hscat=scatterm(splitapply(@mean,data.lat(id),G),splitapply(@mean,data.lon(id),G),[],splitapply(@mean,res.rn(id),G),'filled','MarkerEdgeColor','k');
    %     end
    %     uistack(hcoast,'top')
    
    title(datestr(g.time(mask_fullday(i_t)))); drawnow;
    %     Frame(i_t) = getframe(h);
    %     [imind,cm] = rgb2ind(frame2im(Frame(i_t)),256);
    %     if i_t == 1
    %         imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf);
    %     else
    %         imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append');
    %     end
    delete(hsurf); delete(hscat);
end

v=VideoWriter([filename '.avi']);
v.FrameRate = 4;
v.Quality = 75;
open(v);
writeVideo(v, Frame);
close(v);



%% 3. Reassemble



% real_denst = real_A(:,:,g.day_id(sim_nightly)) + real_I;%nightly_cmt

% if any(real_denst(:)<0)
%     disp(['WARNING: Remove ' num2str(sum(real_denst(:)<0)) ' value(s) because <0'])
%     real_denst(real_denst<0)=nan;
% end

real_dens = nan(size(real_I),'single');

for i_real = 1:nb_real
    real_denst = real_A(:,:,g.day_id,i_real+500) + real_I(:,:,:,i_real);
    real_dens(:,:,:,i_real) = interp1(trans.denst_axis,trans.dens_axis,real_denst,'pchip');
end

%% Save
% save('D:\Guests\rafnuss_at_gmail_com\tmp_BMM2018\Density_simulationMap_reassemble_5','real_dens','-v7.3')
% save('D:\Guests\rafnuss_at_gmail_com\tmp_BMM2018\Density_simulationMap_reassemble_stopover','real_dens','-v7.3') %nightly_cmt
% load('D:\Guests\rafnuss_at_gmail_com\tmp_BMM2018\Density_simulationMap_reassemble.mat')


%% Figure

i_real=1;
dens = real_dens(:,:,:,i_real);

h=figure(2);
filename='data/Density_simulationMap_reassemble';
mask_fullday = find(~reshape(all(all(isnan(dens),1),2),g.nt,1));
worldmap([min(g.lat) max(g.lat)], [min(g.lon) max(g.lon)]);
Frame(numel(mask_fullday)-1) = struct('cdata',[],'colormap',[]);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow('worldrivers.shp','Color', 'blue')
set(gcf,'color','w'); %hcoast=plotm(coastlat, coastlon,'k');
hsurf=[]; hscat=[];
%caxis([0 100]); c=colorbar;c.Label.String = 'Bird density [bird/km^2]';
caxis([0 3]); c=colorbar;c.Label.String = 'Bird density [Log bird/km^2]';

for i_t = 1:numel(mask_fullday)
    
    hsurf=surfm(g.lat2D,g.lon2D,log10(dens(:,:,mask_fullday(i_t))));
    
    %     gtit = datenum(g.time(mask_fullday(i_t)+([-1 0 1])));
    %     id = find(mean(gtit(1:2))<data.time & mean(gtit(2:3))>data.time);
    %     if sum(id)>0
    %         G = findgroups(data.dateradar(id));
    %         hscat=scatterm(splitapply(@mean,data.lat(id),G),splitapply(@mean,data.lon(id),G),[],splitapply(@mean,log(data.dens(id)),G),'filled','MarkerEdgeColor','k');
    %     end
    %     uistack(hcoast,'top')
    
    title(datestr(g.time(mask_fullday(i_t)))); drawnow;
    %     Frame(i_t) = getframe(h);
    %     [imind,cm] = rgb2ind(frame2im(Frame(i_t)),256);
    %     if i_t == 1
    %         imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.1);
    %     else
    %         imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append','DelayTime',0.1);
    %     end
    delete(hsurf);% delete(hscat);
end

v=VideoWriter([filename '.avi']);
v.FrameRate = 4;
v.Quality = 75;
open(v);
writeVideo(v, Frame);
close(v);
