%% Map Estimation of bird mean flight speed
%% 1. Load data

clear all; load('data/dc_corr.mat'); load coastlines; addpath('functions/'); load('data/FlightDir_modelInf.mat')
%% 2. Initialization
%%
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
%% 
% Mask water

g.mask_water = inpolygon(g.lat2D,g.lon2D,coastlat,coastlon);
%% 
% Mask distance Radar

gDdist_sf = pdist2([[dc.lat]' [dc.lon]'], [g.lat2D(:) g.lon2D(:)], @lldistkm);
g.mask_distrad = reshape(min(gDdist_sf)<500,size(g.lat2D));

g.latlonmask = (g.mask_water &  g.mask_distrad);
g.nlm = sum(g.latlonmask(:)==1);
%% 
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
%% 
% Score time of night

g.scoret = (datenum(g.time3D)-mean(cat(4,g.sunrise,g.sunset),4)) ./ (g.sunrise-g.sunset)*2;
g.scoret(g.scoret<-1.1)=-1.1;
g.scoret(g.scoret>1.1)=1.1;
%% 
% Mask of the day

g.mask_day =  g.scoret>-1 & g.scoret<1;
%% 
% Group g.time by night

[~,g.dateradar] = ismember(round(datenum(reshape(g.sunrise(1,1,:),numel(g.time),1))),g.atime);
g.dateradar(g.dateradar==0)=1;
%% 
% Clear the masking variable

clear sunrs_b sunrs_e tim_d lat_d lon_d F 
%% 6.3 Kriging
%%
Gneiting = @(dist,time,range_dist,range_time,delta,gamma,beta) 1./( (time./range_time).^(2.*delta) +1 ) .* exp(-(dist./range_dist).^(2.*gamma)./((time./range_time).^(2.*delta) +1).^(beta.*gamma) );

gtimenum=datenum(g.time);
dir_cov_parm=dir.cov.parm;
g_mask_day=reshape(g.mask_day,g.nl,[]);
g_dateradar=g.dateradar;
ndc = numel(dc);

gDdist_sf = pdist2([[dc.lat]' [dc.lon]'], [g.lat2D(g.latlonmask) g.lon2D(g.latlonmask)], @lldistkm);
g_mask_day = g_mask_day(g.latlonmask(:),:)';

wradius=2;
nb_neigh=100;

gdir_ddn_est=cell(g.nat,1);
gdir_ddn_sig=gdir_ddn_est;
neighday=cell(g.nat,1);
sim=neighday;
dir_cov_C=cell(g.nat,1);
dir_ddn=cell(g.nat,1);
for i_d=1:g.nat
    neighday{i_d}  = find(ismember( data.dateradar, (i_d-1)*ndc+(1:ndc)));
    dir_cov_C{i_d} = dir.cov.C(neighday{i_d},neighday{i_d});
    dir_ddn{i_d} =  dir.ddn(neighday{i_d});
    sim{i_d} = find(g_dateradar==i_d);
    gdir_ddn_est{i_d} = nan(g.nlm,numel(sim{i_d}));
    gdir_ddn_sig{i_d} = gdir_ddn_est{i_d};
end

maxNumCompThreads('automatic');
parpool(maxNumCompThreads());

parfor i_d=1:g.nat
   % maxNumCompThreads(1);
    for i_l=1:g.nlm
        tmpD = gDdist_sf(data.i_r(neighday{i_d}),i_l);
        neighloc = find(tmpD<dir_cov_parm(3)*wradius);
        tmpD=tmpD(neighloc);
        simloc = find(g_mask_day(sim{i_d},i_l));
        for i_sd=1:numel(simloc)
            tmpT = abs(data.time(neighday{i_d}(neighloc))-gtimenum(sim{i_d}(simloc(i_sd))));
            tmpidT = tmpT<dir_cov_parm(4)*wradius;
            neighhour = neighloc(tmpidT);
            if isempty(neighhour)
                gdir_ddn_est{i_d}(i_l,simloc(i_sd)) = 0;
                gdir_ddn_sig{i_d}(i_l,simloc(i_sd)) = sqrt(dir_cov_parm(1) + dir_cov_parm(2));
            else
                Cab = dir_cov_parm(2).*Gneiting(tmpD(tmpidT), tmpT(tmpidT), dir_cov_parm(3), dir_cov_parm(4), dir_cov_parm(5), dir_cov_parm(6), dir_cov_parm(7));
                [Cab,id]=maxk(Cab,nb_neigh);
                lambda = dir_cov_C{i_d}(neighhour(id),neighhour(id))  \  Cab;
                gdir_ddn_est{i_d}(i_l,simloc(i_sd)) = lambda' * dir_ddn{i_d}(neighhour(id));
                gdir_ddn_sig{i_d}(i_l,simloc(i_sd)) = sqrt(dir_cov_parm(1) + dir_cov_parm(2) - lambda' * Cab);
            end
        end
    end
end

gdir.ddn_est=nan(g.nlat,g.nlon,g.nt);
gdir.ddn_sig=gdir.ddn_est;

[latmask, lonmask]=ind2sub([g.nlat g.nlon],latlonmask);
for i_d=1:g.nat
    id = sub2ind([g.nlat g.nlon g.nt], repmat(latmask,numel(sim{i_d}),1), repmat(lonmask,numel(sim{i_d}),1), repelem(sim{i_d},g.nlm));
    gdir.ddn_est(id) = gdir_ddn_est{i_d}(:);
    gdir.ddn_sig(id) = gdir_ddn_sig{i_d}(:);
end
%% 
% Back-transform

gdir.dd_est = gdir.ddn_est+mean(data.dd);
gdir.dd_sig = gdir.ddn_sig;
%% 
% Figure
%%
fig=figure(2);  
filename='data/FlightDir_estimationMap_reassemble';
mask_fullday = find(~reshape(all(all(isnan(gdir.dd_est),1),2),g.nt,1));
Frame(numel(mask_fullday)-1) = struct('cdata',[],'colormap',[]); 
worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]); 
hcoast=plotm(coastlat, coastlon,'k'); caxis([0 360]); hsurf=[];hscat=[];
c=colorbar;c.Label.String = 'Bird Flight Direction [deg]';
for i_t = 1:numel(mask_fullday)
    
    hsurf=surfm(g.lat2D,g.lon2D,gdir.dd_est(:,:,mask_fullday(i_t)));
    
    gtit = datenum(g.time(mask_fullday(i_t)+([-1 0 1])));
    id = find(mean(gtit(1:2))<data.time & mean(gtit(2:3))>data.time);
    if sum(id)>0
        G = findgroups(data.dateradar(id));
        hscat=scatterm(splitapply(@mean,data.lat(id),G),splitapply(@mean,data.lon(id),G),[],splitapply(@mean,data.dd(id),G),'filled','MarkerEdgeColor','k');
    end
    uistack(hcoast,'top')
    
    title(datestr(g.time(mask_fullday(i_t)))); drawnow;

    Frame(i_t) = getframe(fig);
    [imind,cm] = rgb2ind(frame2im(Frame(i_t)),256); 
    if i_t == 1
        imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append');
    end
    delete(hsurf); delete(hscat);
end

v=VideoWriter([filename '.avi']);
v.FrameRate = 4;  % Default 30
v.Quality = 75;
open(v);
writeVideo(v, Frame);
close(v);
%% 
% Save
%%
save('data/FlightDir_estimationMap','g','gampli','gdir','-v7.3')
% load('data/FlightDir_estimationMap')