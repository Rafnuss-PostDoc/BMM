%% Simulation Map
% Script to generate simulated map of bird migration

%% Clear and load data
clear all; load coastlines; addpath('functions/');
load('data/dc_corr');
load('data/FlightDir_modelInf')
load('data/FlightDir_estimationMap');




%% Generation of the path
wradius=2;
neighg_nb=100;
neighr_nb=50;
Gneiting = @(dist,time,range_dist,range_time,delta,gamma,beta) 1./( (time./range_time).^(2.*delta) +1 ) .* exp(-(...
    dist./range_dist).^(2.*gamma)./((time./range_time).^(2.*delta) +1).^(beta.*gamma) );


mask_total = g.mask_day & repmat(g.latlonmask,1,1,g.nt);

ndc = numel(dc);


% 1. Creation of the grid and path
Pathll=cell(1,g.nat);
pathll=cell(1,g.nat);
ndt=cell(1,g.nat);
tic  %3min not-parr
for i_d=1:g.nat
    
    ndt{i_d}=sum(g.dateradar==i_d);
    [LAT,LON,DT]= ndgrid(1:g.nlat,1:g.nlon,1:ndt{i_d});

    Path=Inf*ones(g.nlat,g.nlon,ndt{i_d});
    Path(mask_total(:,:,g.dateradar==i_d))=nan;
    path = nan(sum(isnan(Path(:))),1);
    
    slat = 1:ceil(log(g.nlat+1)/log(2));
    slon = 1:ceil(log(g.nlon+1)/log(2));
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
    
    Pathll{i_d}=reshape(Path(repmat(g.latlonmask,1,1,ndt{i_d})),g.nlm,ndt{i_d});
    [~,pathll{i_d}] = sort(Pathll{i_d}(:));
    pathll{i_d}=pathll{i_d}(~isinf(Pathll{i_d}(pathll{i_d})));
end
toc

% Distance Matrix common to all day
Ddist_gg = squareform(pdist([g.lat2D(g.latlonmask) g.lon2D(g.latlonmask)], @lldistkm));
Ddist_gr = pdist2([[dc.lat]' [dc.lon]'], [g.lat2D(g.latlonmask) g.lon2D(g.latlonmask)], @lldistkm);


%% Computaiton of the weight neigh and S
dir_cov_parm = dir.cov.parm;
NEIGHG=cell(g.nat,1);
NEIGHR=cell(g.nat,1);
LAMBDA=cell(g.nat,1);
S=cell(g.nat,1);

for i_d=1:g.nat
    
    % Sub-selection of the day radar and grid
    neighday{i_d}  = find(ismember( data.dateradar, (i_d-1)*ndc+(1:ndc)));
    dir_cov_C = dir.cov.C(neighday{i_d},neighday{i_d});
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
    
    % 5 Loop of scale for multi-grid path
    tic
    parfor i_pt = 1:numel(pathll{i_d})
        
        % 1. Neighbor from grid
        % Find in the index of the path such that the spatio-temporal distance
        % is within wradius. Then, only select the already simulated value
        % (path value less than currently simulated)
        neighg=find(bsxfun(@and,Ddist_gg(:,LL_i(i_pt))<dir_cov_parm(3)*wradius  , Dtime_gg(:,TT_i(i_pt))'<...
            dir_cov_parm(4)*wradius));
        neighg = neighg(Pathll{i_d}(neighg)<i_pt);
        
        
        Cgp = dir_cov_parm(2).*Gneiting(Ddist_gg(LL_i(Pathll{i_d}(neighg)),LL_i(i_pt)), Dtime_gg(TT_i(Pathll{i_d}(neighg)),...
            TT_i(i_pt)), dir_cov_parm(3), dir_cov_parm(4), dir_cov_parm(5), dir_cov_parm(6), dir_cov_parm(7));
        [Cgp,tmp]=maxk(Cgp,neighg_nb);
        neighg=neighg(tmp);
        
        Cgg = dir_cov_parm(2).*Gneiting(Ddist_gg(LL_i(Pathll{i_d}(neighg)),LL_i(Pathll{i_d}(neighg))), Dtime_gg(TT_i(Pathll{i_d}...
            (neighg)),TT_i(Pathll{i_d}(neighg))), dir_cov_parm(3), dir_cov_parm(4), dir_cov_parm(5), dir_cov_parm(6), dir_cov_parm(7));

        % 2. Find the radar neigh
        
        neighr=find( Ddist_gr(data.i_r(neighday{i_d}),LL_i(i_pt))<dir_cov_parm(3)*wradius & Dtime_gr(:,TT_i(i_pt))<dir_cov_parm(4)*wradius);
        Crp = dir_cov_parm(2).*Gneiting(Ddist_gr(data.i_r(neighday{i_d}(neighr)),LL_i(i_pt)), Dtime_gr(neighr,TT_i(i_pt)), dir_cov_parm(3), ...
            dir_cov_parm(4), dir_cov_parm(5), dir_cov_parm(6), dir_cov_parm(7));
        [Crp,tmp]=maxk(Crp,neighr_nb);
        neighr=neighr(tmp);
        
        Crr = dir_cov_C(neighr,neighr);
        Crg = dir_cov_parm(2).*Gneiting(Ddist_gr(data.i_r(neighday{i_d}(neighr)),LL_i(Pathll{i_d}(neighg))), Dtime_gr(neighr,TT_i(Pathll{i_d}(neighg))), ...
            dir_cov_parm(3), dir_cov_parm(4), dir_cov_parm(5), dir_cov_parm(6), dir_cov_parm(7));
        
        l =  [Cgg Crg'; Crg Crr] \  [Cgp; Crp];
        
        NEIGHG_tmp(:,i_pt) = [neighg ; nan(neighg_nb-numel(neighg),1)];
        NEIGHR_tmp(:,i_pt) = [neighr ; nan(neighr_nb-numel(neighr),1)];
        LAMBDA_tmp(:,i_pt) = [l ; nan(neighr_nb+neighg_nb-numel(l),1)]';
        
        S_tmp(i_pt) = sqrt(dir_cov_parm(1) + dir_cov_parm(2) - l'*[Cgp; Crp]);
        
    end
    
    NEIGHG{i_d} = NEIGHG_tmp;
    NEIGHR{i_d} = NEIGHR_tmp;
    LAMBDA{i_d} = LAMBDA_tmp;
    S{i_d} = S_tmp;
    
    dlmwrite(['data/sim_dir_NEIGHG_' num2str(i_d)],NEIGHG{i_d})
    dlmwrite(['data/sim_dir_NEIGHR_' num2str(i_d)],NEIGHR{i_d})
    dlmwrite(['data/sim_dir_LAMBDA_' num2str(i_d)],LAMBDA{i_d})
    
    toc
    
end

save('data/FlightDir_simulationMap_dir','NEIGHG','NEIGHR','LAMBDA','pathll','S','neighday','sim','-v7.3')


%% 2.2 Residu: Generation of realizations
%load('data/Density_simulationMap_residu')

[latmask, lonmask]=ind2sub([g.nlat g.nlon],g.latlonmask);

nb_real = 1; 
real_dirn = nan(g.nlat,g.nlon,g.nt,nb_real);


for i_real=1:nb_real
    real_dir_ll=nan(g.nlm,g.nt);
    rng('shuffle');
    U=randn(g.nlm,g.nt);
    
    for i_d=1:g.nat
        
        tmp = nan(g.nlm,numel(sim{i_d}));
    
        for i_pt = 1:numel(pathll{i_d})
            ng = ~isnan(NEIGHG{i_d}(:,i_pt));
            nr = ~isnan(NEIGHR{i_d}(:,i_pt));
            nl = ~isnan(LAMBDA{i_d}(:,i_pt));
            tmp(pathll{i_d}(i_pt)) = LAMBDA{i_d}(nl,i_pt)'*[tmp(NEIGHG{i_d}(ng,i_pt)) ; dir.ddn(neighday{i_d}(NEIGHR{i_d}(nr,i_pt)))] + U(i_pt)*S{i_d}(i_pt);
            %assert(~isnan(Resd(pathll{i_d}(i_pt))))
            %assert(isreal(Resd(pathll{i_d}(i_pt))))
        end
        
        real_dir_ll(:,sim{i_d}) = tmp;
    end
    tmp2 = nan(g.nlat,g.nlon,g.nt);
    tmp2(repmat(g.latlonmask,1,1,g.nt))=real_dir_ll;
    real_dirn(:,:,:,i_real) = tmp2;
end

real_dir = real_dirn+mean(data.dd);
real_dir = mod(real_dir,360);

%save('data/FlightDir_simulationMap_reassemble','real_dir')
%load('data/FlightDir_simulationMap_reassemble','real_dir')

i_real=1;

h=figure(2);  
worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]);  
filename='data/FlightDir_simulationMap_reassemble';
mask_fullday = find(~reshape(all(all(isnan(real_dir(:,:,:,i_real)),1),2),g.nt,1));
Frame(numel(mask_fullday)-1) = struct('cdata',[],'colormap',[]); 
hcoast=plotm(coastlat, coastlon,'k'); caxis([0 365]); hscat=[];
for i_t = 1:numel(mask_fullday)

    hsurf=surfm(g.lat2D,g.lon2D,real_dir(:,:,mask_fullday(i_t),i_real));

    gtit = datenum(g.time(mask_fullday(i_t)+([-1 0 1])));
    id = find(mean(gtit(1:2))<data.time & mean(gtit(2:3))>data.time);
    if sum(id)>0
        G = findgroups(data.dateradar(id));
        hscat=scatterm(splitapply(@mean,data.lat(id),G),splitapply(@mean,data.lon(id),G),[],splitapply(@mean,data.dd(id),G),'filled','MarkerEdgeColor','k');
    end
    uistack(hcoast,'top')
    
    title(datestr(g.time(mask_fullday(i_t)))); drawnow;
    Frame(i_t) = getframe(h);
    [imind,cm] = rgb2ind(frame2im(Frame(i_t)),256); 
    if i_t == 1
        imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append');
    end
    delete(hsurf); delete(hscat);
end

v=VideoWriter([filename '.avi']);
v.FrameRate = 4;
v.Quality = 75;
open(v);
writeVideo(v, Frame);
close(v);
