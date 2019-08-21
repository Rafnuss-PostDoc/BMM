%% Map Estimation
% See ... for the model building and parameter inferance.

%% 1. Load data
clear all; 
load('../1-Cleaning/data/dc_corr.mat'); 
load('../1-Cleaning/data/sunrisesunset_grid.mat')
load('../1-Cleaning/data/Rain_grid.mat','TCRW');
load('../2-Inference/data/Density_inference.mat'); 
addpath('./functions/'); 

%% 2. Initialization of the grid
%
g.lat=43:.2:68;
g.lon=-5:.2:30;
g.time = (start_date:1/24/4:end_date);
g.atime=unique(round(datenum(g.time)));

g.nlat=numel(g.lat);
g.nlon=numel(g.lon);
g.nl=numel(g.lat)*numel(g.lon);
g.nt=numel(g.time);
g.nat=numel(g.atime);

[g.lat3D,g.lon3D,g.time3D]=ndgrid(g.lat,g.lon,g.time);
[g.lat2D, g.lon2D] = ndgrid(g.lat,g.lon);

ndc = numel(dc);
dclat = [dc.lat]';
dclon = [dc.lon]';

%% 3. Mask
% Mask water

g.mask_water = inpolygon(g.lat2D,g.lon2D,coastlat,coastlon);

% figure; hold on; 
% worldmap([min(dclat) max(dclat)], [min(dclon) max(dclon)]);  
% plotm(g.lat2D(g.mask_water),g.lon2D(g.mask_water),'xr')
% plotm(g.lat2D(~g.mask_water),g.lon2D(~g.mask_water),'xb')
% plotm(coastlat, coastlon,'k')
% legend({'Land','Water'})

% Mask distance Radar
Ddist = pdist2([dclat dclon], [g.lat2D(:) g.lon2D(:)], @lldistkm);
g.mask_distrad = reshape(min(Ddist)<250,size(g.lat2D));

g.latlonmask = (g.mask_water &  g.mask_distrad);
g.nlm = sum(g.latlonmask(:));

% figure; hold on; 
% worldmap([min(dclat) max(dclat)], [min(dclon) max(dclon)]);  
% plotm(g.lat2D(g.mask_distrad),g.lon2D(g.mask_distrad),'xr')
% plotm(g.lat2D(~g.mask_distrad),g.lon2D(~g.mask_distrad),'xb')
% plotm(dclat', dclon','ok')
% plotm(coastlat, coastlon,'k')

% Sunrise and sunset
g.sunrise = nan(g.nlat,g.nlon,g.nt);
g.sunset = g.sunrise;
for i_t=1:numel(tim_d)-1
    id = g.time>tim_d(i_t) & g.time<=tim_d(i_t+1);
    
    F=griddedInterpolant({lat_d,lon_d},sunrs_e(:,:,i_t));
    g.sunset(:,:,id)=repmat(F(g.lat2D,g.lon2D),1,1,sum(id));
    
    F=griddedInterpolant({lat_d,lon_d},sunrs_b(:,:,i_t+1));
    g.sunrise(:,:,id)=repmat(F(g.lat2D,g.lon2D),1,1,sum(id));
end

% figure('position',[0 0 1000 400]); hold on;
% plot(g.time(:),reshape(mean(mean(g.sunset,1),2),1,2017))
% plot(g.time(:),reshape(mean(mean(g.sunrise,1),2),1,2017))
% legend('sunset','sunrise')

% Score time of night
g.NNT = (datenum(g.time3D)-mean(cat(4,g.sunrise,g.sunset),4)) ./ (g.sunrise-g.sunset)*2;
%g.NNT(g.NNT<-1.1)=-1.1;
%g.NNT(g.NNT>1.1)=1.1;

% filename='data/';
% Frame(numel(g.time)-1) = struct('cdata',[],'colormap',[]); 
% fig=figure(2);
% worldmap([min(dclat) max(dclat)], [min(dclon) max(dclon)]); 
% hcoast=plotm(coastlat, coastlon,'k'); hsurf=[];
% for i_t = 1:numel(g.time)-1
%     delete(hsurf);
%     hsurf=surfm(g.lat2D,g.lon2D,g.NNT(:,:,i_t));
%     uistack(hcoast,'top')
%     title(datestr(g.time(i_t))); caxis([-1 1]); drawnow;
%     Frame(i_t) = getframe(fig);
%     [imind,cm] = rgb2ind(frame2im(Frame(i_t)),256); 
%     if i_t == 1
%         imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append');
%     end
% end
% v=VideoWriter([filename '.avi']);
% v.FrameRate = 4;
% v.Quality = 75;
% open(v);
% writeVideo(v, Frame);
% close(v);

% figure;
% plot(g.time3D(:), g.NNT(:), '.k' )

% Mask of the day
g.mask_day =  g.NNT>-1 & g.NNT<1;

% Group g.time by night
[~,g.dateradar] = ismember(round(reshape(g.sunrise(1,1,:),g.nt,1)),g.atime);
% g.dateradar(g.dateradar==0)=1;

% Mask Rain
rain = [];
for i_d=1:ndc
    id=dc(i_d).scoret<1&dc(i_d).scoret>-1;
    rain=[rain; repmat(i_d,sum(id),1) dc(i_d).denss(id), dc(i_d).tcrw(id) nansum(dc(i_d).DBZH(id,:),2)];
end

% figure;
% subplot(3,1,1);scatter(rain(:,2),rain(:,4),'.'); xlabel('density');ylabel('DBZH')
% subplot(3,1,2);scatter(rain(:,2),rain(:,3),'.');xlabel('density');ylabel('tcrw')
% subplot(3,1,3);scatter(rain(:,3),rain(:,4),'.');xlabel('tcrw');ylabel('DBZH')
% 
% figure; hold on; plot(rain(:,2)); yyaxis right; plot(rain(:,3))
% 
% % histogram of rain
% histogram(log10(rain(:,3)))
% 
% histogram(log10(rain(:,2)))
% 
% % Compare histogram of rain vs no rain
% id=rain(:,3)<10^-8;
% figure; hold on;
% histogram((rain(id,2)))
% histogram((rain(~id,2)))
% legend('no Rain','Rain')
% 
% VAL=0:0.001:0.1;
% r=nan(size(VAL));
% r2=nan(size(VAL));
% r3=nan(size(VAL));
% for i_val=1:numel(VAL)
%     id=rain(:,3)>VAL(i_val);
%     r(i_val)=mean(rain(id,2)==0);
%     r2(i_val)=mean(rain(~id,2)==0);
%     r3(i_val)=sum(id);
% end
% figure; hold on; plot(VAL,r);plot(VAL,r2); 
%  yyaxis right; plot(VAL,r3); 
% xlabel('threashold'); ylabel('percent of 0 value');
% legend('in removed data','in kept data')

% Threashold of 0.03 (0.02-0.05)
g.mask_rain_thr=0.03;
g.mask_rain=TCRW<g.mask_rain_thr;

% Clear the masking variable
clear sunrs_b sunrs_e lat_d lon_d F tim_d dc


%% 6.2 Amplitude Kriging

% Covariance matrix of weather radars
Dtime=squareform(pdist(ampli.T(~ampli.isnan)));
Ddist=squareform(pdist([dclat(ampli.R(~ampli.isnan)) dclon(ampli.R(~ampli.isnan))],@lldistkm));
Caa = ampli.cov.C(Ddist,Dtime);

% Cross-covariance of the radar data
Ddist = pdist2([dclat dclon], [g.lat2D(g.latlonmask) g.lon2D(g.latlonmask)], @lldistkm);
Ddist_sf = Ddist(ampli.R(~ampli.isnan),repmat((1:g.nlm)',g.nat,1));
Dtime_sf = pdist2(ampli.T(~ampli.isnan), repelem(g.atime',g.nlm,1));

Cab = ampli.cov.C(Ddist_sf,Dtime_sf);

% Compute the Kriging Weight
Lambda = inv(Caa)*Cab;

% Kriging estimation
gampli.A_est = nan(g.nlat,g.nlon,g.nat);
gampli.A_est(repmat(g.latlonmask,1,1,g.nat)) = (Lambda' * ampli.An) + repmat(g.lat2D(g.latlonmask),g.nat,1)*ampli.w(2)+ampli.w(1);

% Kriging Variance
gampli.A_sig = nan(g.nlat,g.nlon,g.nat);
gampli.A_sig(repmat(g.latlonmask,1,1,g.nat)) = sqrt(sum(ampli.cov.parm(1:3)) - sum(Lambda.*Cab));



%% Figures
figure('position',[0 0 1000 600]); 
for i_t = 1:g.nat-1
    subplot(4,6,i_t); hold on; 
    h=worldmap([min(dclat) max(dclat)], [min(dclon) max(dclon)]);  
    setm(h,'frame','on','grid','off'); set(findall(h,'Tag','MLabel'),'visible','off'); set(findall(h,'Tag','PLabel'),'visible','off')
    geoshow('landareas.shp', 'FaceColor', [215 215 215]./255); geoshow('worldrivers.shp','Color', 'blue')
    surfm(g.lat2D,g.lon2D,gampli.A_est(:,:,i_t)); caxis([.6 1.9])
    scatterm(dclat,dclon,[],ampli.A(:,i_t),'filled','MarkerEdgeColor','k')
end

% filename='data/Density_estimationMap_amplitude_A';
filename='data/Density_estimationMap_amplitude_An';
Frame(numel(g.atime)-1) = struct('cdata',[],'colormap',[]); 
fig=figure(2); 
h=worldmap([min(dclat) max(dclat)], [min(dclon) max(dclon)]);  
setm(h,'frame','on','grid','off'); set(findall(h,'Tag','MLabel'),'visible','off'); set(findall(h,'Tag','PLabel'),'visible','off')
geoshow('landareas.shp', 'FaceColor', [215 215 215]./255); geoshow('worldrivers.shp','Color', 'blue')
set(gcf,'color','w'); hsurf=[];hscat=[];
caxis([-.5 .5]); 
%caxis([0 80]);  c=colorbar;c.Label.String = 'Bird density [bird/km^2]';c.Label.FontSize=14;
for i_t = 1:numel(g.atime)-1
    delete(hsurf); delete(hscat)
    
    hsurf=surfm(g.lat2D,g.lon2D,gampli.A_est(:,:,i_t)-g.lat2D*ampli.w(2)-ampli.w(1));
    hscat=scatterm(dclat,dclon,[],ampli.A(:,i_t)-dclat(ampli.r)*ampli.w(2)-ampli.w(1),'filled','MarkerEdgeColor','k');
%     hsurf=surfm(g.lat2D,g.lon2D,gampli.A_est(:,:,i_t).^(1/pow_a));
%     hscat=scatterm(dclat,dclon,[],ampli.A(:,i_t).^(1/pow_a),'filled','MarkerEdgeColor','k');
    
    title(datestr(g.atime(i_t)),'FontSize',14); drawnow;
    Frame(i_t) = getframe(fig);
    [imind,cm] = rgb2ind(frame2im(Frame(i_t)),256); 
    if i_t == 1
        imwrite(imind,cm,[filename '.gif'],'gif','Loopcount',inf);
    else
        imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append');
    end
end


%% 6.3 Residu Kriging
% Covariance matrix of weather radar data 
Dtime = sq80uareform(pdist(data.time));
Ddist_sm = squareform(pdist([dclat dclon], @lldistkm));
Ddist = Ddist_sm(data.radar,data.radar);
tic;Caa = res.cov.C(Ddist,Dtime);toc

Dtime_sm = single(pdist2(data.time, datenum(g.time')));
Ddist_sm = single(pdist2([dclat dclon], [g.lat2D(:) g.lon2D(:)], @lldistkm));

gres_Rn_est = nan(size(g.mask_day));
gres_Rn_sig = nan(size(g.mask_day));

for i_t=1:g.nat
    mask_r = data.day==i_t;

    mask_g = g.mask_day & repmat(reshape(g.dateradar==i_t,1,1,[]),g.nlat,g.nlon,1) & repmat(g.latlonmask,1,1,g.nt);
    
    [tmp_i,tmp_j,tmp_k] = ind2sub(size(g.time3D),find(mask_g));
    
    Dtime_sf = Dtime_sm(mask_r,tmp_k);
    Ddist_sf = Ddist_sm(data.radar(mask_r), sub2ind(size(g.lat2D),tmp_i,tmp_j));
    
    tic;Cab = res.cov.C(Ddist_sf,Dtime_sf);toc

    Lambda = inv(Caa(mask_r,mask_r))*Cab;

    gres_Rn_est(mask_g) = Lambda' * res.Rn(mask_r);
    gres_Rn_sig(mask_g) = sqrt(sum(res.cov.parm(1:3)) - sum(Lambda.*Cab)); 
    
end

gres.R_est = gres_Rn_est .* sqrt(polyval(flipud(res.b),g.NNT)) + polyval(flipud(res.a),g.NNT);
gres.R_sig = sqrt( gres_Rn_sig.^2 + polyval(flipud(res.b),g.NNT) );


%% Figure
% filename='data/Density_estimationMap_residu_R';
filename='data/Density_estimationMap_residu_Rn';
mask_fullday = find(~reshape(all(all(isnan(gres_Rn_est),1),2),g.nt,1));
Frame(numel(mask_fullday)-1) = struct('cdata',[],'colormap',[]);
fig=figure;
h=worldmap([min(dclat) max(dclat)], [min(dclon) max(dclon)]); 
setm(h,'grid','off'); set(findall(h,'Tag','MLabel'),'visible','off'); set(findall(h,'Tag','PLabel'),'visible','off')
geoshow('landareas.shp', 'FaceColor', [215 215 215]./255); geoshow('worldrivers.shp','Color', 'blue')
set(gcf,'color','w'); hsurf=[];hscat=[];
caxis([-5 5]);

for i_t = 2:numel(mask_fullday)-1
    delete(hsurf); delete(hscat);
    hsurf=surfm(g.lat2D,g.lon2D,gres_Rn_est(:,:,mask_fullday(i_t)));
    gtit = datenum(g.time(mask_fullday(i_t)+([-1 0 1])));
    id = find(mean(gtit(1:2))<data.time & mean(gtit(2:3))>data.time);
    if sum(id)>0
        [G,ID] = findgroups(data.radar(id));
        hscat=scatterm(dclat(ID),dclon(ID),[],splitapply(@mean,res.Rn(id),G),'filled','MarkerEdgeColor','k');
    end
    title(datestr(g.time(mask_fullday(i_t))),'FontSize',14); drawnow
    Frame(i_t) = getframe(fig);
    [imind,cm] = rgb2ind(frame2im(Frame(i_t)),256); 
    if i_t == 2
        imwrite(imind,cm,[filename '.gif'],'gif', 'DelayTime',0.1,'Loopcount',inf);
    else
        imwrite(imind,cm,[filename '.gif'],'gif', 'DelayTime',0.1,'WriteMode','append');
    end
end
% v=VideoWriter([filename '.avi']);
% v.FrameRate = 4;
% v.Quality = 75;
% open(v);
% writeVideo(v, Frame);
% close(v);

%% 6.4 Reassemble
%

gd.denst_est = gampli.A_est(:,:,g.dateradar) + gres.R_est;
gd.denst_sig = sqrt(gampli.A_sig(:,:,g.dateradar).^2 + gres.R_sig.^2);

if any(gd.denst_est(:)<0)
    disp(['WARNING: Remove ' num2str(sum(gd.denst_est(:)<0)) ' value(s) because <0'])
    gd.denst_est(gd.denst_est<0)=nan;
end

gd.dens_est = (gd.denst_est).^(1/pow_a);

gd.dens_q10 = ( gd.denst_est+norminv(.1).*gd.denst_sig ).^(1/pow_a);
gd.dens_q90 = ( gd.denst_est+norminv(.9).*gd.denst_sig ).^(1/pow_a);


%% Save
save('data/Density_estimationMap','g','gampli','gres','gd','-v7.3');
% load('data/Density_estimationMap')


%% Figure
fig=figure(2);  
filename='data/Density_estimationMap_reassamble_Z';
% mask_fullday = 1553:1609;
filename='data/Density_estimationMap_reassamble';
mask_fullday = find(~reshape(all(all(isnan(gd.dens_est),1),2),g.nt,1));
Frame(numel(mask_fullday)-1) = struct('cdata',[],'colormap',[]); 
h=worldmap([min(dclat) max(dclat)], [min(dclon) max(dclon)]); 
setm(h,'grid','off'); set(findall(h,'Tag','MLabel'),'visible','off'); set(findall(h,'Tag','PLabel'),'visible','off')
geoshow('landareas.shp', 'FaceColor', [215 215 215]./255); geoshow('worldrivers.shp','Color', 'blue')
set(gcf,'color','w'); %hcoast=plotm(coastlat, coastlon,'k');
hsurf=[]; hscat=[];
caxis([0 200]); c=colorbar;c.Label.String = 'Bird density [bird/km^2]';c.Label.FontSize=14;
for i_t = 2:numel(mask_fullday)
    delete(hsurf); delete(hscat);
    hsurf=surfm(g.lat2D,g.lon2D,gd.dens_est(:,:,mask_fullday(i_t)));
    gtit = datenum(g.time(mask_fullday(i_t)+([-1 0 1])));
    id = find(mean(gtit(1:2))<data.time & mean(gtit(2:3))>data.time);
    if sum(id)>0
        [G,ID] = findgroups(data.radar(id));
        hscat=scatterm(dclat(ID),dclon(ID),[],splitapply(@mean,data.dens(id),G),'filled','MarkerEdgeColor','k');
    end

    title(datestr(g.time(mask_fullday(i_t))),'FontSize',14); drawnow;

    Frame(i_t) = getframe(fig);
    [imind,cm] = rgb2ind(frame2im(Frame(i_t)),256); 
    if i_t == 2
        imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append','DelayTime',0.1);
    end
    delete(hsurf);% delete(hscat);
end


