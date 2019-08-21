%% Simulation Map
% Script to generate simulated map of bird migration

%% Clear and load data
clear all; 
load('../1-Cleaning/data/dc_corr.mat'); 
load('../2-Inference/data/Density_inference.mat'); 
load('../4-Estimation/data/Density_estimationMap','g');
load coastlines;
addpath('./functions/'); 

dclat = [dc.lat]' ;
dclon = [dc.lon]';

%% 1. Amplitude

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


% FIGURE
% rLAT=repmat(dclat,1,g.nat); rLAT=rLAT(~ampli.isnan);
% rLON=repmat(dclon,1,g.nat); rLON=rLON(~ampli.isnan);
% rT = repmat(g.atime,numel(dc),1); rT=rT(~ampli.isnan);
% LAT_tmp = LAT(repmat(g.latlonmask,1,1,g.nat));
% LON_tmp = LON(repmat(g.latlonmask,1,1,g.nat));
% AT_tmp = AT(repmat(g.latlonmask,1,1,g.nat));
% 
% figure(30);clf; hold on; view(3);
% % imagesc(g.lat,g.lon,g.latlonmask')
% % scatter3(rLAT(:),rLON(:),rT(:),'filled');
% for i_pt=1:numel(pathll)
%     scatter3(g.lat(LAT_tmp(pathll(i_pt))),g.lon(LON_tmp(pathll(i_pt))),g.atime(AT_tmp(pathll(i_pt))),[],i_pt,'filled')
% %     if mod(i_pt,10)==0, keyboard; end
% %     keyboard
% end

% 1.2 Ditance Matrix
Ddist_gg = squareform(pdist([g.lat2D(g.latlonmask) g.lon2D(g.latlonmask)], @lldistkm));
Ddist_gr = pdist2(repmat([dclat dclon], g.nat,1), [g.lat2D(g.latlonmask) g.lon2D(g.latlonmask)], @lldistkm);
Ddist_gr = Ddist_gr(~ampli.isnan,:);
Dtime_gg = squareform(pdist(g.atime'));
Dtime_gr = pdist2(repelem(g.atime',numel(dc),1), g.atime');
Dtime_gr = Dtime_gr(~ampli.isnan,:);
Dtime_rr=squareform(pdist(ampli.T(~ampli.isnan)));
Ddist_rr=squareform(pdist([dclat(ampli.R(~ampli.isnan)) dclon(ampli.R(~ampli.isnan))],@lldistkm));
C_rr = ampli.cov.C(Ddist_rr,Dtime_rr);

% 1.3. Initizialization of the kriging weights and variance error,Computation of neighbours, lambda weight and S
neighg_nb=100;
neighr_nb=50;
NEIGHG = nan(neighg_nb,numel(path),'single');
NEIGHR = nan(neighr_nb,numel(path),'single');
LAMBDA = nan(neighg_nb+neighr_nb,numel(path),'single');
S = nan(numel(path),1,'single');
tmp = repmat((1:g.nlm)',1,g.nat);
LL_i=tmp(pathll);
tmp = repmat(1:numel(g.atime),g.nlm,1);
TT_i=tmp(pathll);

wradius=3;


parfor i_pt = 1:numel(path)
    
    % 1. Neighbor from grid
    % Find in the index of the path such that the spatio-temporal distance
    % is within wradius. Then, only select the already simulated value
    % (path value less than currently simulated)
    neighg=find(Pathll<i_pt & bsxfun(@and,Ddist_gg(:,LL_i(i_pt))<ampli.cov.parm(4)*wradius  , Dtime_gg(:,TT_i(i_pt))'<ampli.cov.parm(5)*wradius));    
    Cgp = ampli.cov.C(Ddist_gg(LL_i(Pathll(neighg)),LL_i(i_pt)), Dtime_gg(TT_i(Pathll(neighg)),TT_i(i_pt)));
    [Cgp,tmp]=maxk(Cgp,neighg_nb);
    neighg=neighg(tmp);
    Cgg = ampli.cov.C(Ddist_gg(LL_i(Pathll(neighg)),LL_i(Pathll(neighg))), Dtime_gg(TT_i(Pathll(neighg)),TT_i(Pathll(neighg))));

    
    % 2. Find the radar neigh
    neighr=find(Ddist_gr(:,LL_i(i_pt))<ampli.cov.parm(4)*wradius & Dtime_gr(:,TT_i(i_pt))<ampli.cov.parm(5)*wradius);
    Crp = ampli.cov.C(Ddist_gr(neighr,LL_i(i_pt)), Dtime_gr(neighr,TT_i(i_pt)));
    [Crp,tmp]=maxk(Crp,neighr_nb);
    neighr=neighr(tmp);
    Crr = C_rr(neighr,neighr);
    
    Crg = ampli.cov.C(Ddist_gr(neighr,LL_i(Pathll(neighg))), Dtime_gr(neighr,TT_i(Pathll(neighg))));
    
    l =  [Cgg Crg'; Crg Crr] \ [Cgp; Crp];
    
    NEIGHG(:,i_pt) = [neighg ; nan(neighg_nb-numel(neighg),1)];
    NEIGHR(:,i_pt) = [neighr ; nan(neighr_nb-numel(neighr),1)];
    LAMBDA(:,i_pt) = [l ; nan(neighr_nb+neighg_nb-numel(l),1)]';
    
    S(i_pt) = sqrt(sum(ampli.cov.parm(1:3)) - l'*[Cgp; Crp]);
end


% save('data/Density_simulationMap_amplitude','NEIGHG','NEIGHR','LAMBDA','pathll','S')
% load('data/Density_simulationMap_amplitude')


%% 2. Amplitude: Simulation of realization

% load('data/Density_simulationMap_amplitude')

nb_real = 100;
real_An_ll = nan(g.nlm,g.nat,nb_real);
ampli_An = ampli.An;

for i_real=1:nb_real
    real_ampli_tmp=nan(g.nlm,g.nat);
    rng('shuffle');
    U=randn(g.nlm,g.nat);
    
    for i_pt = 1:numel(pathll)
        ng = ~isnan(NEIGHG(:,i_pt));
        nr = ~isnan(NEIGHR(:,i_pt));
        nl = ~isnan(LAMBDA(:,i_pt));
        real_ampli_tmp(pathll(i_pt)) = LAMBDA(nl,i_pt)'*[real_ampli_tmp(NEIGHG(ng,i_pt)) ; ampli_An(NEIGHR(nr,i_pt))] + U(i_pt)*S(i_pt);
        assert(~isnan(real_ampli_tmp(pathll(i_pt))))
        assert(isreal(real_ampli_tmp(pathll(i_pt))))
    end
    real_An_ll(:,:,i_real) = real_ampli_tmp;
end

real_A = nan(g.nlat,g.nlon,g.nat,nb_real);
real_A(repmat(g.latlonmask,1,1,g.nat,nb_real)) = real_An_ll + repmat(g.lat2D(g.latlonmask),1,g.nat,nb_real)*ampli.w(2) + ampli.w(1);

%save('real_A_100','real_A')

figure('position',[0 0 1000 600]);
for i_t = 1:g.nat-1
    subplot(4,6,i_t); hold on; worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]);
    surfm(g.lat2D,g.lon2D,real_A(:,:,i_t))
    scatterm([dc.lat],[dc.lon],[],ampli.A(:,i_t)','filled','MarkerEdgeColor','k')
    plotm(coastlat, coastlon,'k')
end

% fig=figure(2);  
% filename='data/Density_simulationMap_amplitude';
% Frame(g.nat-1) = struct('cdata',[],'colormap',[]); 
% worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]); 
% hcoast=plotm(coastlat, coastlon,'k'); caxis([min(ampli.A(:)) max(ampli.A(:))]);
% for i_t = 1:g.nat-1
%      
%     hsurf=surfm(g.lat2D,g.lon2D,real_A(:,:,i_t));
%     hscat=scatterm([dc.lat],[dc.lon],[],ampli.A(:,i_t),'filled','MarkerEdgeColor','k');
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
mask_total = g.mask_day & repmat(g.latlonmask,1,1,g.nt);
ndc = numel(dc);

Pathll=cell(1,g.nat);
pathll=cell(1,g.nat);
ndt=cell(1,g.nat);

slat = 1:ceil(log(g.nlat+1)/log(2));
slon = 1:ceil(log(g.nlon+1)/log(2));
    
tic  %3min not-parr because of ismember (lin 220)
for i_d=1:g.nat
    ndt{i_d}=sum(g.dateradar==i_d);
    [LAT,LON,DT]= ndgrid(1:g.nlat,1:g.nlon,1:ndt{i_d});

    Path=Inf*ones(g.nlat,g.nlon,ndt{i_d});
    Path(mask_total(:,:,g.dateradar==i_d))=nan;
    path = nan(sum(isnan(Path(:))),1);
    
    sat = 1:ceil(log(ndt{i_d}+1)/log(2));
    sn = max([numel(slat), numel(slon) numel(sat)]);
    nb = nan(sn,1);
    start = zeros(sn+1,1);
    ds = 2.^(sn-1:-1:0);
    for i_scale = 1:sn
        [LAT_s,LON_s,DT_s] = ndgrid(1:ds(i_scale):g.nlat,1:ds(i_scale):g.nlon,1:ds(i_scale):ndt{i_d}); % matrix coordinate
        id = find(isnan(Path(:)));
        id = id(ismember([LAT(id) LON(id) DT(id)], [LAT_s(:) LON_s(:) DT_s(:)], 'rows'));
        nb(i_scale) = numel(id);
        start(i_scale+1) = start(i_scale)+nb(i_scale);
        path( start(i_scale)+(1:nb(i_scale)) ) = id(randperm(nb(i_scale)));
        Path(path( start(i_scale)+(1:nb(i_scale)) )) = start(i_scale)+(1:nb(i_scale));
    end
    
    Pathll{i_d}=reshape(Path(repmat(g.latlonmask,1,1,ndt{i_d})), g.nlm, ndt{i_d});
    [~,pathll{i_d}] = sort(Pathll{i_d}(:));
    pathll{i_d}=pathll{i_d}(~isinf(Pathll{i_d}(pathll{i_d})));
    
    Pathll{i_d}= int64(Pathll{i_d});
    pathll{i_d}= int64(pathll{i_d});
end
toc

% Distance Matrix common to all day
Ddist_gg = squareform(pdist([g.lat2D(g.latlonmask) g.lon2D(g.latlonmask)], @lldistkm));
Ddist_gr = pdist2([dclat dclon], [g.lat2D(g.latlonmask) g.lon2D(g.latlonmask)], @lldistkm);

% Covariance Matrix of hard data
Dtime = squareform(pdist(data.time));
Ddist_sm = squareform(pdist([dclat dclon], @lldistkm));
Ddist = Ddist_sm(data.radar,data.radar);
tic;Crrn = res.cov.C(Ddist,Dtime);toc

%% 
% 2.2 Residu: Computaiton of the weight neigh and S
res_cov_parm = res.cov.parm;
NEIGHG=cell(g.nat,1);
NEIGHR=cell(g.nat,1);
LAMBDA=cell(g.nat,1);
S=cell(g.nat,1);
neighday=cell(g.nat,1);
res_cov_C=cell(g.nat,1);
sim=cell(g.nat,1);

wradius=4;
neighg_nb=100;
neighr_nb=50;


for i_d=1:g.nat
    
    % Sub-selection of the day radar and grid
    neighday{i_d}  = find(data.day==i_d);
    res_cov_C_i_d = Crrn(neighday{i_d},neighday{i_d});
    sim{i_d} = find(g.dateradar==i_d);
   
    % Distance matrix for the day
    Dtime_gg = squareform(pdist(datenum(g.time(sim{i_d}))'));
    Dtime_gr = pdist2(data.time(neighday{i_d}), datenum(g.time(sim{i_d}))');
    
    % 4. Initizialization of the kriging weights and variance error
    NEIGHG_tmp = nan(neighg_nb,numel(pathll{i_d}));
    NEIGHR_tmp = nan(neighr_nb,numel(pathll{i_d}));
    LAMBDA_tmp = nan(neighg_nb+neighr_nb,numel(pathll{i_d}));
    S_tmp = nan(numel(pathll{i_d}),1);
    tmp = repmat((1:g.nlm)',1,ndt{i_d});
    LL_i=tmp(pathll{i_d});
    tmp = repmat(1:ndt{i_d},g.nlm,1);
    TT_i=tmp(pathll{i_d});
    
    Pathll_i_d = Pathll{i_d};
    
    % 5 Loop of scale for multi-grid path
    tic
    parfor i_pt = 1:numel(pathll{i_d})
        
        % 1. Neighbor from grid
        % Find in the index of the path such that the spatio-temporal distance
        % is within wradius. Then, only select the already simulated value
        % (path value less than currently simulated)
        neighg = find( Pathll_i_d<i_pt & bsxfun(@and,Ddist_gg(:,LL_i(i_pt))<res_cov_parm(4)*wradius  , Dtime_gg(:,TT_i(i_pt))'<res_cov_parm(5)*wradius));
        Cgp = res.cov.C(Ddist_gg(LL_i(Pathll_i_d(neighg)),LL_i(i_pt)), Dtime_gg(TT_i(Pathll_i_d(neighg)), TT_i(i_pt)));
        [Cgp,tmp]=maxk(Cgp,neighg_nb);
        neighg=neighg(tmp);
        
        Cgg = res.cov.C(Ddist_gg(LL_i(Pathll_i_d(neighg)),LL_i(Pathll_i_d(neighg))), Dtime_gg(TT_i(Pathll_i_d(neighg)),TT_i(Pathll_i_d(neighg))));

        % 2. Find the radar neigh
        neighr=find( Ddist_gr(data.radar(neighday{i_d}),LL_i(i_pt))<res_cov_parm(4)*wradius & Dtime_gr(:,TT_i(i_pt))<res_cov_parm(5)*wradius);
        Crp = res.cov.C(Ddist_gr(data.radar(neighday{i_d}(neighr)),LL_i(i_pt)), Dtime_gr(neighr,TT_i(i_pt)));
        [Crp,tmp]=maxk(Crp,neighr_nb);
        neighr=neighr(tmp);
        Crr = res_cov_C_i_d(neighr,neighr);
        
        % 3. cross covariance
        Crg = res.cov.C(Ddist_gr(data.radar(neighday{i_d}(neighr)),LL_i(Pathll_i_d(neighg))), Dtime_gr(neighr,TT_i(Pathll_i_d(neighg))));
        
        l =  [Cgg Crg'; Crg Crr] \  [Cgp; Crp];
        
        NEIGHG_tmp(:,i_pt) = [neighg ; nan(neighg_nb-numel(neighg),1)];
        NEIGHR_tmp(:,i_pt) = [neighr ; nan(neighr_nb-numel(neighr),1)];
        LAMBDA_tmp(:,i_pt) = [l ; nan(neighr_nb+neighg_nb-numel(l),1)]';
        
        S_tmp(i_pt) = sqrt(sum(res_cov_parm(1:3)) - l'*[Cgp; Crp]);

    end
    
    NEIGHG{i_d} = single(NEIGHG_tmp);
    NEIGHR{i_d} = single(NEIGHR_tmp);
    LAMBDA{i_d} = single(LAMBDA_tmp);
    S{i_d} = single(S_tmp);
    
%     dlmwrite(['data/sim_res_NEIGHG_' num2str(i_d)],NEIGHG{i_d})
%     dlmwrite(['data/sim_res_NEIGHR_' num2str(i_d)],NEIGHR{i_d})
%     dlmwrite(['data/sim_res_LAMBDA_' num2str(i_d)],LAMBDA{i_d})
    
    toc
    
end


% save('data/Density_simulationMap_residu','NEIGHG','NEIGHR','LAMBDA','pathll','S','neighday','sim','-v7.3')
% load('data/Density_simulationMap_residu')

%% 2.2 Residu: Generation of realizations
%load('data/Density_simulationMap_residu')

[latmask, lonmask]=ind2sub([g.nlat g.nlon],find(g.latlonmask));

res_Rn = res.Rn;
nb_real = 100; 
real_resn = nan(g.nlat,g.nlon,g.nt,nb_real,'single');


for i_real=1:nb_real
   
    real_res_ll=nan(g.nlm,g.nt);
    rng('shuffle');
    U=randn(g.nlm,g.nt);
    
    for i_d=1:g.nat
        
        tmp = nan(g.nlm,numel(sim{i_d}));
    
        for i_pt = 1:numel(pathll{i_d})
            ng = ~isnan(NEIGHG{i_d}(:,i_pt));
            nr = ~isnan(NEIGHR{i_d}(:,i_pt));
            nl = ~isnan(LAMBDA{i_d}(:,i_pt));
            tmp(pathll{i_d}(i_pt)) = LAMBDA{i_d}(nl,i_pt)'*[tmp(NEIGHG{i_d}(ng,i_pt)) ; res_Rn(neighday{i_d}(NEIGHR{i_d}(nr,i_pt)))] + U(i_pt)*S{i_d}(i_pt);
            %assert(~isnan(Resd(pathll{i_d}(i_pt))))
            %assert(isreal(Resd(pathll{i_d}(i_pt))))
        end
        
        real_res_ll(:,sim{i_d}) = tmp;
    end
    tmp2 = nan(g.nlat,g.nlon,g.nt);
    tmp2(repmat(g.latlonmask,1,1,g.nt))= real_res_ll;
    real_resn(:,:,:,i_real) = tmp2;

end

real_res = real_resn .* sqrt(polyval(flipud(res.b),repmat(g.NNT,1,1,1,nb_real))) + polyval(flipud(res.a),repmat(g.NNT,1,1,1,nb_real));


i_real=1;

h=figure(2);  
worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]);  
filename='data/Density_simulationMap_residu';
mask_fullday = find(~reshape(all(all(isnan(real_res(:,:,:,i_real)),1),2),g.nt,1));
Frame(numel(mask_fullday)-1) = struct('cdata',[],'colormap',[]); 
hcoast=plotm(coastlat, coastlon,'k'); caxis([-3 3]); hscat=[];
for i_t = 2:numel(mask_fullday)

    hsurf=surfm(g.lat2D,g.lon2D,real_resn(:,:,mask_fullday(i_t),i_real));

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


real_denst = real_A(:,:,g.dateradar) + real_res;

if any(real_denst(:)<0)
    disp(['WARNING: Remove ' num2str(sum(real_denst(:)<0)) ' value(s) because <0'])
    real_denst(real_denst<0)=nan;
end

real_dens = (real_denst).^(1/pow_a);




% save('data/Density_simulationMap_reassemble','real_dens','-v7.3')
% load('data/Density_simulationMap_reassemble.mat')

real_dens_ll=reshape(real_dens(repmat(g.latlonmask,1,1,g.nt,nb_real)),g.nlm,g.nt,nb_real);

% save('data/Density_simulationMap_reassemble_ll','real_dens','-v7.3')
% load('data/Density_simulationMap_reassemble_ll.mat')

%% Figure

i_real=1;
dens = nan(g.nlat,g.nlon,g.nt);
dens(repmat(g.latlonmask,1,1,g.nt)) = real_dens_ll(:,:,i_real);

h=figure(2);  
filename='data/Density_simulationMap_reassemble';
mask_fullday = find(~reshape(all(all(isnan(dens),1),2),g.nt,1));
worldmap([min(g.lat) max(g.lat)], [min(g.lon) max(g.lon)]);  
Frame(numel(mask_fullday)-1) = struct('cdata',[],'colormap',[]); 
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow('worldrivers.shp','Color', 'blue')
set(gcf,'color','w'); %hcoast=plotm(coastlat, coastlon,'k');
hsurf=[]; hscat=[];caxis([0 3]); c=colorbar;c.Label.String = 'Bird density [Log bird/km^2]';

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
    Frame(i_t) = getframe(h);
    [imind,cm] = rgb2ind(frame2im(Frame(i_t)),256); 
    if i_t == 1
        imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append','DelayTime',0.1);
    end
    delete(hsurf);% delete(hscat);
end

v=VideoWriter([filename '.avi']);
v.FrameRate = 4;
v.Quality = 75;
open(v);
writeVideo(v, Frame);
close(v);


%% Check sim vs est

load('data/Density_estimationMap')

% err sim
% ok if transformed, slight asymetry if not. 
err = reshape(gd.dens_est(repmat(g.latlonmask,1,1,g.nt)),sum(g.latlonmask(:)),g.nt).^pow_a - mean(real_dens_ll.^pow_a,3);
histogram(err(:))
err = reshape(gd.dens_est(repmat(g.latlonmask,1,1,g.nt)),sum(g.latlonmask(:)),g.nt) - mean(real_dens_ll,3);
histogram(err(:))

real_dens_ll_s=sort(real_dens_ll,3);

err = reshape(gd.dens_q10(repmat(g.latlonmask,1,1,g.nt)),sum(g.latlonmask(:)),g.nt)-real_dens_ll_s(:,:,10);
histogram(err(:))
err = reshape(gd.dens_q90(repmat(g.latlonmask,1,1,g.nt)),sum(g.latlonmask(:)),g.nt)-real_dens_ll_s(:,:,90);
histogram(err(:))

