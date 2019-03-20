%% Simulation Map
% Script to generate simulated map of bird migration

%% Clear and load data
clear all; load coastlines; addpath('../functions/');
load('data/dc_corr');
load('data/Flight_modelInf')
load('data/Flight_estimationMap');
load('data/Density_estimationMap.mat','g')




%% Compute neighborhood, weight and variance 

% Generation of the path
mask_total = g.mask_day & repmat(g.latlonmask,1,1,g.nt);
ndc = numel(dc);

Pathll=cell(1,g.nat);
pathll=cell(1,g.nat);
ndt=cell(1,g.nat);
tic  %3min not-parr because of ismember (lin 220)
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
    
    Pathll{i_d}= int64(Pathll{i_d});
    pathll{i_d}= int64(pathll{i_d});
end
toc

% Distance Matrix common to all day
Ddist_gg = squareform(pdist([g.lat2D(g.latlonmask) g.lon2D(g.latlonmask)], @lldistkm));
Ddist_gr = pdist2([[dc.lat]' [dc.lon]'], [g.lat2D(g.latlonmask) g.lon2D(g.latlonmask)], @lldistkm);


% 2.2 Residu: Computaiton of the weight neigh and S
uv_cov_parm_u=uv.cov.parm_u;
uv_cov_parm_v=uv.cov.parm_v;
NEIGHG=cell(g.nat,1);
NEIGHR=cell(g.nat,1);
LAMBDA_u=cell(g.nat,1);
LAMBDA_v=cell(g.nat,1);
S_u=cell(g.nat,1);
S_v=cell(g.nat,1);

dist_thr_g = 500;
time_thr_g = 0.4;
dist_thr_r = 2000;
time_thr_r = 2;
neighg_nb=100;
neighr_nb=50;
Gneiting = @(dist,time,range_dist,range_time,delta,gamma,beta) 1./( (time./range_time).^(2.*delta) +1 ) .* exp(-(...
    dist./range_dist).^(2.*gamma)./((time./range_time).^(2.*delta) +1).^(beta.*gamma) );

Cf_u = @(dist,time) uv_cov_parm_u(2).*Gneiting( dist, time, uv_cov_parm_u(3), uv_cov_parm_u(4), uv_cov_parm_u(5), uv_cov_parm_u(6), uv_cov_parm_u(7) );
Cf_v = @(dist,time) uv_cov_parm_v(2).*Gneiting( dist, time, uv_cov_parm_v(3), uv_cov_parm_v(4), uv_cov_parm_v(5), uv_cov_parm_v(6), uv_cov_parm_v(7) );


for i_d=1:g.nat
    % Sub-selection of the day radar and grid
    neighday{i_d}  = find(ismember( data.dateradar, (i_d-1)*ndc+(1:ndc)));
    uv_cov_C_u_i_d = uv.cov.C_u(neighday{i_d},neighday{i_d});
    uv_cov_C_v_i_d = uv.cov.C_v(neighday{i_d},neighday{i_d});
    uv_utrans{i_d} =  uv.utrans(neighday{i_d});
    uv_vtrans{i_d} =  uv.vtrans(neighday{i_d});
    sim{i_d} = find(g.dateradar==i_d);
    
    % Distance matrix for the day
    Dtime_gg = squareform(pdist(datenum(g.time(sim{i_d}))'));
    Dtime_gr = pdist2(data.time(neighday{i_d}), datenum(g.time(sim{i_d}))');
    
    % 4. Initizialization of the kriging weights and variance error
    NEIGHG_tmp = nan(neighg_nb,numel(pathll{i_d}));
    NEIGHR_tmp = nan(neighr_nb,numel(pathll{i_d}));
    LAMBDA_u_tmp = nan(neighg_nb+neighr_nb,numel(pathll{i_d}));
    LAMBDA_v_tmp = nan(neighg_nb+neighr_nb,numel(pathll{i_d}));
    S_u_tmp = nan(numel(pathll{i_d}),1);
    S_v_tmp = nan(numel(pathll{i_d}),1);
    tmp_u = repmat((1:g.nlm)',1,ndt{i_d});
    LL_i=tmp_u(pathll{i_d});
    tmp_u = repmat(1:ndt{i_d},g.nlm,1);
    TT_i=tmp_u(pathll{i_d});
    Pathll_i_d = Pathll{i_d};
    neighday_i_d = neighday{i_d};
    
    % 5 Loop of scale for multi-grid path
    tic
    parfor i_pt = 1:numel(pathll{i_d})
        
        % 1. Neighbor from grid
        % Find in the index of the path such that the spatio-temporal distance
        % is within wradius. Then, only select the already simulated value
        % (path value less than currently simulated)
        neighg = find(Pathll_i_d<i_pt & bsxfun(@and, Ddist_gg(:,LL_i(i_pt))<dist_thr_g , Dtime_gg(:,TT_i(i_pt))'<time_thr_g ));
        
        Cgp_u = Cf_u( Ddist_gg(LL_i(Pathll_i_d(neighg)),LL_i(i_pt)), Dtime_gg(TT_i(Pathll_i_d(neighg)),TT_i(i_pt)) );
        [Cgp_u,tmp_u]=maxk(Cgp_u,neighg_nb);
        neighg=neighg(tmp_u);
        
        Cgg_u = uv_cov_parm_u(1)*eye(numel(neighg)) + Cf_u(Ddist_gg(LL_i(Pathll_i_d(neighg)),LL_i(Pathll_i_d(neighg))), Dtime_gg(TT_i(Pathll_i_d(neighg)),TT_i(Pathll_i_d(neighg))));
        Cgp_v = Cf_v(Ddist_gg(LL_i(Pathll_i_d(neighg)),LL_i(i_pt)), Dtime_gg(TT_i(Pathll_i_d(neighg)),TT_i(i_pt)));
        Cgg_v = uv_cov_parm_v(1)*eye(numel(neighg)) + Cf_v(Ddist_gg(LL_i(Pathll_i_d(neighg)),LL_i(Pathll_i_d(neighg))), Dtime_gg(TT_i(Pathll_i_d(neighg)),TT_i(Pathll_i_d(neighg))));
        
        
        % 2. Find the radar neigh
        neighr=find( Ddist_gr(data.i_r(neighday_i_d),LL_i(i_pt))<dist_thr_r & Dtime_gr(:,TT_i(i_pt))<time_thr_r);
        Crp_u = Cf_u(Ddist_gr(data.i_r(neighday_i_d(neighr)),LL_i(i_pt)), Dtime_gr(neighr,TT_i(i_pt)));
        [Crp_u,tmp_u]=maxk(Crp_u,neighr_nb);
        neighr=neighr(tmp_u);
        Crp_v = Cf_v(Ddist_gr(data.i_r(neighday_i_d(neighr)),LL_i(i_pt)), Dtime_gr(neighr,TT_i(i_pt)));
        
        Crr_u = uv_cov_C_u_i_d(neighr,neighr);
        Crr_v = uv_cov_C_v_i_d(neighr,neighr);
        Crg_u = Cf_u(Ddist_gr(data.i_r(neighday_i_d(neighr)),LL_i(Pathll_i_d(neighg))), Dtime_gr(neighr,TT_i(Pathll_i_d(neighg))));
        Crg_v = Cf_v(Ddist_gr(data.i_r(neighday_i_d(neighr)),LL_i(Pathll_i_d(neighg))), Dtime_gr(neighr,TT_i(Pathll_i_d(neighg))));
        
        l_u =  [Cgg_u Crg_u'; Crg_u Crr_u] \  [Cgp_u; Crp_u];
        l_v =  [Cgg_v Crg_v'; Crg_v Crr_v] \  [Cgp_v; Crp_v];
        
        NEIGHG_tmp(:,i_pt) = [neighg ; nan(neighg_nb-numel(neighg),1)];
        NEIGHR_tmp(:,i_pt) = [neighr ; nan(neighr_nb-numel(neighr),1)];
        
        LAMBDA_u_tmp(:,i_pt) = [l_u ; nan(neighr_nb+neighg_nb-numel(l_u),1)]';
        LAMBDA_v_tmp(:,i_pt) = [l_v ; nan(neighr_nb+neighg_nb-numel(l_v),1)]';
        S_u_tmp(i_pt) = sqrt(uv_cov_parm_u(1) + uv_cov_parm_u(2) - l_u'*[Cgp_u; Crp_u]);
        S_v_tmp(i_pt) = sqrt(uv_cov_parm_v(1) + uv_cov_parm_v(2) - l_v'*[Cgp_v; Crp_v]);
        assert(isreal(S_v_tmp(i_pt)))
    end
    
    NEIGHG{i_d} = single(NEIGHG_tmp);
    NEIGHR{i_d} = single(NEIGHR_tmp);
    LAMBDA_u{i_d} = single(LAMBDA_u_tmp);
    LAMBDA_v{i_d} = single(LAMBDA_v_tmp);
    S_u{i_d} = single(S_u_tmp);
    S_v{i_d} = single(S_v_tmp);
    
%     dlmwrite(['data/sim_flight_NEIGHG_' num2str(i_d)],NEIGHG{i_d})
%     dlmwrite(['data/sim_flight_NEIGHR_' num2str(i_d)],NEIGHR{i_d})
%     dlmwrite(['data/sim_flight_LAMBDA_u_' num2str(i_d)],LAMBDA_u{i_d})
%     dlmwrite(['data/sim_flight_LAMBDA_v_' num2str(i_d)],LAMBDA_v{i_d})
%     dlmwrite(['data/sim_flight_S_u_' num2str(i_d)],S_u{i_d})
%     dlmwrite(['data/sim_flight_S_v_' num2str(i_d)],S_v{i_d})    
    toc
    
end


for i_d=1:g.nat
    NEIGHG{i_d} = single(NEIGHG{i_d} );
    NEIGHR{i_d} = single(NEIGHR{i_d});
    LAMBDA_u{i_d} = single(LAMBDA_u{i_d});
    LAMBDA_v{i_d} = single(LAMBDA_v{i_d});
    S_u{i_d} = single(S_u{i_d});
    S_v{i_d} = single(S_v{i_d});
end


% save('data/Flight_simulationMap','NEIGHG','NEIGHR','LAMBDA_u','LAMBDA_v','pathll','S_u','S_v','neighday','sim','-v7.3')
% load('data/Flight_simulationMap')

%% 2.2 Generation of realizations
%load('data/Density_simulationMap_residu')

[latmask, lonmask]=ind2sub([g.nlat g.nlon],find(g.latlonmask));

nb_real = 100; 
real_un_ll = nan(g.nlm,g.nt,nb_real);
real_vn_ll = nan(g.nlm,g.nt,nb_real);


for i_real=1:nb_real
    rng('shuffle');
    Uu=randn(g.nlm,g.nt);
    Uv=randn(g.nlm,g.nt);
    
    for i_d=1:g.nat
        
        tmp_u = nan(g.nlm,numel(sim{i_d}));
        tmp_v = nan(g.nlm,numel(sim{i_d}));
    
        for i_pt = 1:numel(pathll{i_d})
            ng = ~isnan(NEIGHG{i_d}(:,i_pt));
            nr = ~isnan(NEIGHR{i_d}(:,i_pt));
            nl = ~isnan(LAMBDA_u{i_d}(:,i_pt));
            tmp_u(pathll{i_d}(i_pt)) = LAMBDA_u{i_d}(nl,i_pt)'*[tmp_u(NEIGHG{i_d}(ng,i_pt)) ; uv.utrans(neighday{i_d}(NEIGHR{i_d}(nr,i_pt)))] + Uu(i_pt)*S_u{i_d}(i_pt);
            tmp_v(pathll{i_d}(i_pt)) = LAMBDA_v{i_d}(nl,i_pt)'*[tmp_v(NEIGHG{i_d}(ng,i_pt)) ; uv.vtrans(neighday{i_d}(NEIGHR{i_d}(nr,i_pt)))] + Uu(i_pt)*S_v{i_d}(i_pt);
            %assert(~isnan(Resd(pathll{i_d}(i_pt))))
            %assert(isreal(Resd(pathll{i_d}(i_pt))))
        end
        
        real_un_ll(:,sim{i_d},i_real) = tmp_u;
        real_vn_ll(:,sim{i_d},i_real) = tmp_v;
        
    end
end

real_u_ll = single(real_un_ll *uv.trans.std(1) + uv.trans.mean(1));
real_v_ll = single(real_vn_ll *uv.trans.std(2) + uv.trans.mean(2));


save('data/Flight_simulationMap_real_ll','real_u_ll','real_v_ll','-v7.3');

i_real=1;

h=figure(2);  
worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]);  
filename='data/Density_simulationMap_residu';
mask_fullday = find(~reshape(all(all(isnan(real_u(:,:,:,i_real)),1),2),g.nt,1));
Frame(numel(mask_fullday)-1) = struct('cdata',[],'colormap',[]); 
hcoast=plotm(coastlat, coastlon,'k'); caxis([-3 3]); hscat=[];
for i_t = 1:numel(mask_fullday)

    hsurf=quiverm(g.lat2D,g.lon2D,real_u(:,:,mask_fullday(i_t),i_real),real_v(:,:,mask_fullday(i_t),i_real));

%     gtit = datenum(g.time(mask_fullday(i_t)+([-1 0 1])));
%     id = find(mean(gtit(1:2))<data.time & mean(gtit(2:3))>data.time);
%     if sum(id)>0
%         G = findgroups(data.dateradar(id));
%         hscat=scatterm(splitapply(@mean,data.lat(id),G),splitapply(@mean,data.lon(id),G),[],splitapply(@mean,res.rn(id),G),'filled','MarkerEdgeColor','k');
%     end
%     uistack(hcoast,'top')
    drawnow
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
