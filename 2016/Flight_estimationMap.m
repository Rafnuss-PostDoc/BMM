%% Map Estimation

%% 1. Load data
%
clear all; load('data/dc_corr.mat'); load coastlines; addpath('../functions/'); load('data/Flight_modelInf.mat')


%% 2. Initialization of the grid
%
g.lat=43:.2:68;
g.lon=-5:.2:30;
g.time = start_date-1:1/24/4:end_date-1;
g.atime=unique(round(datenum(g.time)));

g.nlat=numel(g.lat);
g.nlon=numel(g.lon);
g.nl=numel(g.lat)*numel(g.lon);
g.nt=numel(g.time);
g.nat=numel(g.atime);

[g.lat3D,g.lon3D,g.time3D]=ndgrid(g.lat,g.lon,g.time);
[g.lat2D, g.lon2D] = ndgrid(g.lat,g.lon);

%% 3. Mask
% 
% Mask water
g.mask_water = inpolygon(g.lat2D,g.lon2D,coastlat,coastlon);

figure; hold on; 
worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]);  
plotm(g.lat2D(g.mask_water),g.lon2D(g.mask_water),'xr')
plotm(g.lat2D(~g.mask_water),g.lon2D(~g.mask_water),'xb')
plotm(coastlat, coastlon,'k')
legend({'Land','Water'})
 
% Mask distance Radar
max_dist = 200; % corresponding to 1.2 * ampli.cov.parm(3)
gDdist_sf = pdist2([[dc.lat]' [dc.lon]'], [g.lat2D(:) g.lon2D(:)], @lldistkm);
g.mask_distrad = reshape(min(gDdist_sf)<max_dist,size(g.lat2D));
g.latlonmask = (g.mask_water &  g.mask_distrad);
g.nlm = sum(g.latlonmask(:));

figure; hold on; 
worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]);  
plotm(g.lat2D(g.mask_distrad),g.lon2D(g.mask_distrad),'xr')
plotm(g.lat2D(~g.mask_distrad),g.lon2D(~g.mask_distrad),'xb')
plotm([dc.lat]', [dc.lon]','ok')
plotm(coastlat, coastlon,'k')

% Sunrise and sunset
load('./data/sunrisesunset_grid.mat')
g.sunrise = nan(g.nlat,g.nlon,g.nt);
g.sunset = g.sunrise;
for i_t=2:numel(tim_d)
    id = g.time>tim_d(i_t)-1 & g.time<=tim_d(i_t);
    
    F=griddedInterpolant({lat_d,lon_d},sunrs_e(:,:,i_t-1));
    g.sunset(:,:,id)=repmat(F(g.lat2D,g.lon2D),1,1,sum(id));
    
    F=griddedInterpolant({lat_d,lon_d},sunrs_b(:,:,i_t));
    g.sunrise(:,:,id)=repmat(F(g.lat2D,g.lon2D),1,1,sum(id));
end
figure('position',[0 0 1000 400]); hold on;
plot(g.time(:),reshape(mean(mean(g.sunset,1),2),1,2017))
plot(g.time(:),reshape(mean(mean(g.sunrise,1),2),1,2017))
legend('sunset','sunrise')

% Score time of night
g.scoret = (datenum(g.time3D)-mean(cat(4,g.sunrise,g.sunset),4)) ./ (g.sunrise-g.sunset)*2;
g.scoret(g.scoret<-1.1)=-1.1;
g.scoret(g.scoret>1.1)=1.1;
figure; plot(g.time3D(:), g.scoret(:), '.k' )

% Mask of the day
g.mask_day =  g.scoret>-1 & g.scoret<1;

% Group g.time by night
[~,g.dateradar] = ismember(round(datenum(reshape(g.sunrise(1,1,:),numel(g.time),1))),g.atime);
g.dateradar(g.dateradar==0)=1;

% Mask Rain
% Threashold of 0.03 (0.02-0.05) see Density_estimationMap for calculation
g.mask_rain_thr=0.03;
load('data/Rain_grid.mat','TCRW');
g.mask_rain=TCRW<0.03;

% Clear the masking variable
clear sunrs_b sunrs_e tim_d lat_d lon_d F 


%% 6.3 Residu Kriging
%
Gneiting = @(dist,time,range_dist,range_time,delta,gamma,beta) 1./( (time./range_time).^(2.*delta) +1 ) .* exp(-(dist./range_dist).^(2.*gamma)./((time./range_time).^(2.*delta) +1).^(beta.*gamma) );

gtimenum=datenum(g.time);
uv_cov_parm_u=uv.cov.parm_u;
uv_cov_parm_v=uv.cov.parm_v;
g_mask_day=reshape(g.mask_day,g.nl,[]);
g_dateradar=g.dateradar;
ndc = numel(dc);


gDdist_sf = pdist2([[dc.lat]' [dc.lon]'], [g.lat2D(g.latlonmask) g.lon2D(g.latlonmask)], @lldistkm);
g_mask_day = g_mask_day(g.latlonmask(:),:)';

wradius=2;
neigh_nb=100;

guv_trans_est=cell(g.nat,1);
guv_trans_sig=guv_trans_est;
neighday=cell(g.nat,1);
sim=neighday;
uv_cov_C_u=cell(g.nat,1);
uv_utrans=cell(g.nat,1);
uv_vtrans=cell(g.nat,1);
for i_d=1:g.nat
    neighday{i_d}  = find(ismember( data.dateradar, (i_d-1)*ndc+(1:ndc)));
    uv_cov_C_u{i_d} = uv.cov.C_u(neighday{i_d},neighday{i_d});
    uv_cov_C_v{i_d} = uv.cov.C_v(neighday{i_d},neighday{i_d});
    uv_utrans{i_d} =  uv.utrans(neighday{i_d});
    uv_vtrans{i_d} =  uv.vtrans(neighday{i_d});
    sim{i_d} = find(g_dateradar==i_d);
    guv_utrans_est{i_d} = nan(g.nlm,numel(sim{i_d}));
    guv_utrans_sig{i_d} = guv_utrans_est{i_d};
    guv_vtrans_est{i_d} = guv_utrans_est{i_d};
    guv_vtrans_sig{i_d} = guv_utrans_est{i_d};
end

% maxNumCompThreads('automatic');
% parpool(maxNumCompThreads());

for i_d=1:g.nat
   %maxNumCompThreads(1);
    for i_l=1:g.nlm
        tmpD = gDdist_sf(data.i_r(neighday{i_d}),i_l);
        neighloc = find(tmpD<uv_cov_parm_v(3)*wradius);
        tmpD=tmpD(neighloc);
        simloc = find(g_mask_day(sim{i_d},i_l));
        for i_sd=1:numel(simloc)
            tmpT = abs(data.time(neighday{i_d}(neighloc))-gtimenum(sim{i_d}(simloc(i_sd))));
            tmpidT = tmpT<uv_cov_parm_v(4)*wradius;
            neighhour = neighloc(tmpidT);
            if isempty(neighhour)
                guv_utrans_est{i_d}(i_l,simloc(i_sd)) = 0;
                guv_utrans_sig{i_d}(i_l,simloc(i_sd)) = uv_cov_parm_u(1) + uv_cov_parm_u(2);
                guv_vtrans_est{i_d}(i_l,simloc(i_sd)) = 0;
                guv_vtrans_sig{i_d}(i_l,simloc(i_sd)) = uv_cov_parm_v(1) + uv_cov_parm_v(2);
            else
                Cab_u = uv_cov_parm_u(2).*Gneiting(tmpD(tmpidT), tmpT(tmpidT), uv_cov_parm_u(3), uv_cov_parm_u(4), uv_cov_parm_u(5), uv_cov_parm_u(6), uv_cov_parm_u(7));
                [~,id]=mink(Cab_u,neigh_nb);
                lambda = uv_cov_C_u{i_d}(neighhour(id),neighhour(id))  \  Cab_u(id);
                guv_utrans_est{i_d}(i_l,simloc(i_sd)) = lambda' * uv_utrans{i_d}(neighhour(id));
                guv_utrans_sig{i_d}(i_l,simloc(i_sd)) = sqrt(uv_cov_parm_u(1) + uv_cov_parm_u(2) - lambda' * Cab_u(id));

                Cab_v = uv_cov_parm_v(2).*Gneiting(tmpD(tmpidT), tmpT(tmpidT), uv_cov_parm_v(3), uv_cov_parm_v(4), uv_cov_parm_v(5), uv_cov_parm_v(6), uv_cov_parm_v(7));
                [~,id]=mink(Cab_v,neigh_nb);
                lambda = uv_cov_C_v{i_d}(neighhour(id),neighhour(id))  \  Cab_v(id);
                guv_vtrans_est{i_d}(i_l,simloc(i_sd)) = lambda' * uv_vtrans{i_d}(neighhour(id));
                guv_vtrans_sig{i_d}(i_l,simloc(i_sd)) = sqrt(uv_cov_parm_v(1) + uv_cov_parm_v(2) - lambda' * Cab_v(id));
            end
        end
    end
end

guv.utrans_est=nan(g.nlat,g.nlon,g.nt);
guv.vtrans_est=nan(g.nlat,g.nlon,g.nt);
guv.vtrans_sig=nan(g.nlat,g.nlon,g.nt);
guv.vtrans_sig=nan(g.nlat,g.nlon,g.nt);

[latmask, lonmask]=ind2sub([g.nlat g.nlon],find(g.latlonmask));
for i_d=1:g.nat
    id = sub2ind([g.nlat g.nlon g.nt], repmat(latmask,numel(sim{i_d}),1), repmat(lonmask,numel(sim{i_d}),1), repelem(sim{i_d},g.nlm));
    guv.utrans_est(id) = guv_utrans_est{i_d}(:);
    guv.utrans_sig(id) = guv_utrans_sig{i_d}(:);
    guv.vtrans_est(id) = guv_vtrans_est{i_d}(:);
    guv.vtrans_sig(id) = guv_vtrans_sig{i_d}(:);
end

 
%% Back-transform
%
a = uv.tran.finv([guv.utrans_est guv.vtrans_est],uv);
guv.u_est=a(:,1);
guv.v_est=a(:,2);

a = uv.tran.finv([guv.utrans_est+norminv(.1).*guv.utrans_sig guv.vtrans_est+norminv(.1).*guv.vtrans_sig],uv);
guv.u_p10=a(:,1);
guv.v_p10=a(:,2);

a = uv.tran.finv([guv.utrans_est+norminv(.9).*guv.utrans_sig guv.vtrans_est+norminv(.9).*guv.vtrans_sig],uv);
guv.u_p90=a(:,1);
guv.v_p90=a(:,2);

%% Figure
%
filename='data/Density_estimationMap_residu';
mask_fullday = find(~reshape(all(all(isnan(gres.rn_est),1),2),g.nt,1));
Frame(numel(mask_fullday)-1) = struct('cdata',[],'colormap',[]);

fig=figure(2);
worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]); 
hcoast=plotm(coastlat, coastlon,'k'); hsurf=[]; hscat=[];
for i_t = 1:numel(mask_fullday)
    delete(hsurf); delete(hscat);
    hsurf=surfm(g.lat2D,g.lon2D,gres.r_est(:,:,mask_fullday(i_t)));
    gtit = datenum(g.time(mask_fullday(i_t)+([-1 0 1])));
    id = find(mean(gtit(1:2))<data.time & mean(gtit(2:3))>data.time);
    if sum(id)>0
        G = findgroups(data.dateradar(id));
        hscat=scatterm(splitapply(@mean,data.lat(id),G),splitapply(@mean,data.lon(id),G),[],splitapply(@mean,res.r(id),G),'filled','MarkerEdgeColor','k');
    end
    uistack(hcoast,'top');
    caxis([-.26 .27])
    title(datestr(g.time(mask_fullday(i_t)))); drawnow
    Frame(i_t) = getframe(fig);
    [imind,cm] = rgb2ind(frame2im(Frame(i_t)),256); 
    if i_t == 1
        imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append');
    end
end
v=VideoWriter([filename '.avi']);
v.FrameRate = 4;
v.Quality = 75;
open(v);
writeVideo(v, Frame);
close(v);


% Figure
fig=figure(2);  
filename='data/Density_estimationMap_reassamble_05102016';
mask_fullday = find(~reshape(all(all(isnan(g.dens_est),1),2),g.nt,1));
mask_fullday = 1553:1609;
Frame(numel(mask_fullday)-1) = struct('cdata',[],'colormap',[]); 
worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]);  
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow('worldrivers.shp','Color', 'blue')
set(gcf,'color','w'); %hcoast=plotm(coastlat, coastlon,'k');
 hsurf=[]; hscat=[];
caxis([-2 6]); c=colorbar;c.Label.String = 'Bird density [Log bird/m^2]';
for i_t = 1:numel(mask_fullday)
    
    hsurf=surfm(g.lat2D,g.lon2D,log(g.dens_est(:,:,mask_fullday(i_t))));
    gtit = datenum(g.time(mask_fullday(i_t)+([-1 0 1])));
    id = find(mean(gtit(1:2))<data.time & mean(gtit(2:3))>data.time);
%     if sum(id)>0
%         G = findgroups(data.dateradar(id));
%         hscat=scatterm(splitapply(@mean,data.lat(id),G),splitapply(@mean,data.lon(id),G),[],splitapply(@mean,log(data.dens(id)),G),'filled','MarkerEdgeColor','k');
%     end
    % uistack(hcoast,'top');
    title(datestr(g.time(mask_fullday(i_t)))); drawnow;

    Frame(i_t) = getframe(fig);
    [imind,cm] = rgb2ind(frame2im(Frame(i_t)),256); 
    if i_t == 1
        imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append','DelayTime',0.1);
    end
    delete(hsurf);% delete(hscat);
end

v=VideoWriter([filename '.avi']);
v.FrameRate = 4;  % Default 30
v.Quality = 75;
open(v);
writeVideo(v, Frame);
close(v);

%%  Save
%
save('data/Flight_estimationMap','g','guv','-v7.3')
% load('data/Flight_estimationMap')