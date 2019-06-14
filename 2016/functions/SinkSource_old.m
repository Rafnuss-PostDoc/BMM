%% Sink/source version with padding on the side of 0. Not ideal

% Estimation 
load('data/Density_estimationMap','gd','g')
load('data/Flight_estimationMap','guv')
addpath('./functions/')
load coastlines;

dn = 2;
gext.lat = fillmissing([nan(1,dn) g.lat nan(1,dn)],'linear');
gext.lon = fillmissing([nan(1,dn) g.lon nan(1,dn)],'linear');
gext.nlat = numel(gext.lat);
gext.nlon = numel(gext.lon);
[gext.lon2D, gext.lat2D] = meshgrid(gext.lon,gext.lat);

vlon = zeros(gext.nlat,gext.nlon,g.nt);
vlon(dn+1:end-dn,dn+1:end-dn,:) = guv.u_est; % m/s (+) east, (-) west
vlon(isnan(vlon))=0;

vlat = zeros(gext.nlat,gext.nlon,g.nt);
vlat(dn+1:end-dn,dn+1:end-dn,:) = guv.v_est; % m/s (+) east, (-) west
vlat(isnan(vlat))=0;

c = zeros(gext.nlat,gext.nlon,g.nt);
tmp=gd.dens_est;
tmp(~g.mask_rain)=0;
c(dn+1:end-dn,dn+1:end-dn,:) = gd.dens_est; % m/s (+) east, (-) west
c(isnan(c))=0;

% compute area
dlat2 = lldistkm([gext.lat(1) gext.lon(1)],[gext.lat(2) gext.lon(1)]);
dlon2 = lldistkm([gext.lat2D(:,1) gext.lon2D(:,1)],[gext.lat2D(:,1) gext.lon2D(:,2)]);
area = repmat(dlon2*dlat2,1,gext.nlon,g.nt);
gext.latlonmask = padarray(double(g.latlonmask),[dn dn],2);

% Define the delta distance
dlat = lldistkm([gext.lat(1) 0],[gext.lat(2) 0]) * ones(gext.nlat,gext.nlon)*1000;
dlon = repmat(lldistkm([gext.lat(:) gext.lon(1)*ones(gext.nlat,1)],[gext.lat(:) gext.lon(2)*ones(gext.nlat,1)]), 1, gext.nlon)*1000;
dt=15*60;

% Central Difference with variable \Delta x
% dxdlon = @(x) 1/2 * ((x(2:end-1,3:end,:)-x(2:end-1,2:end-1,:))./repmat(dlon(2:end-1,2:end),1,1,g.nt) + (x(2:end-1,2:end-1,:)-x(2:end-1,1:end-2,:))./repmat(dlon(2:end-1,1:end-1),1,1,g.nt));
% dxdlat = @(x) 1/2 * ((x(3:end,2:end-1,:)-x(2:end-1,2:end-1,:))./repmat(dlat(2:end,2:end-1),1,1,g.nt) + (x(2:end-1,2:end-1,:)-x(1:end-2,2:end-1,:))./repmat(dlat(1:end-1,2:end-1),1,1,g.nt));
dxdlon = @(x) 1/2 * (( [x(:,2:end,:) zeros(size(x,1),1,size(x,3))]-x)./repmat(dlon,1,1,g.nt) + (x-[zeros(size(x,1),1,size(x,3)) x(:,1:end-1,:)])./repmat(dlon,1,1,g.nt));
dxdlat = @(x) 1/2 * (( [x(2:end,:,:); zeros(1,size(x,2),size(x,3))]-x)./repmat(dlat,1,1,g.nt) + (x-[zeros(1,size(x,2),size(x,3)); x(1:end-1,:,:)])./repmat(dlat,1,1,g.nt));


W = nan(gext.nlat,gext.nlon,g.nt);

dcvdlon = dxdlon(c.*vlon);
dcvdlon2D = dcvdlon;
dcvdlon2D(:,2:end-1,:) = 1/4*( dcvdlon(:,1:end-2,:) +2*dcvdlon(:,2:end-1,:) +dcvdlon(:,3:end,:) ); 
dcvdlat = dxdlat(c.*vlat);
dcvdlat2D = dcvdlat;
dcvdlat2D(2:end-1,:,:) = 1/4*( dcvdlat(1:end-2,:,:) +2*dcvdlat(2:end-1,:,:) +dcvdlat(3:end,:,:) ); 
F = dcvdlon2D+dcvdlat2D;
c_tp1 = c(:,:,1:end-1) - 0.5*dt*(F(:,:,2:end)+F(:,:,1:end-1));
W(:,:,1:end-1) = c(:,:,2:end)-c_tp1;

% W(~repmat(gext.latlonmask,1,1,g.nt) & W==0)=nan;
 W(W==0)=nan;

% save('data/SinkSourceSim.mat','W')











% Map of total time average
W_all=nansum(W(:,:,:),3);
W_all(W_all==0)=NaN;

% PROBLEM !!!!
W_all_2 =W_all;
W_all_2(gext.latlonmask==2)=nan;

figure('position',[0 0 1000 600]); 
colormap(brewermap([],'Spectral'))
worldmap([gext.lat(1) gext.lat(end)], [gext.lon(1) gext.lon(end)]);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow('worldrivers.shp','Color', 'blue')
surfm(gext.lat2D,gext.lon2D,W_all); 
caxis([min(W_all(:)) max(W_all(:))])
col=colorbar;col.Label.String = 'Bird Sink/source [bird/km^2]';
plotm(coastlat, coastlon,'k')



% Positiv and negativ flux
W_arrival=W; W_departure=W;
W_arrival(W_arrival>=0)=0;
W_departure(W_departure<=0)=0;




% Maps of daily average
W_day=nan(gext.nlat,gext.nlon,g.nat);
W_day_departure=nan(gext.nlat,gext.nlon,g.nat);
W_day_arrival=nan(gext.nlat,gext.nlon,g.nat);
for i_t=1:g.nat-1
    idt=g.atime(i_t)-0.5<datenum(g.time) & g.atime(i_t)+0.5>datenum(g.time);
    W_day(:,:,i_t)=nansum(W(:,:,idt),3);
    W_day_departure(:,:,i_t)=nansum(W_departure(:,:,idt),3);
    W_day_arrival(:,:,i_t)=nansum(W_arrival(:,:,idt),3);
end
W_day(W_day==0)=nan;
W_day_departure(W_day_departure==0)=nan;
W_day_arrival(W_day_arrival==0)=nan;

figure('position',[0 0 1000 600]); 
colormap(brewermap([],'Spectral'))
for i_t = 1:g.nat-1
    subplot(4,6,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
    surfm(gext.lat2D,gext.lon2D,W_day(:,:,i_t)); caxis([-100 100])
    plotm(coastlat, coastlon,'k')
end

figure('position',[0 0 1000 600]); 
for i_t = 1:g.nat-1
    subplot(4,6,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
    surfm(gext.lat2D,gext.lon2D,W_day_departure(:,:,i_t)); caxis([0 100])
    plotm(coastlat, coastlon,'k')
end

figure('position',[0 0 1000 600]); 
for i_t = 1:g.nat-1
    subplot(4,6,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
    surfm(gext.lat2D,gext.lon2D,W_day_arrival(:,:,i_t)); caxis([-100 0])
    plotm(coastlat, coastlon,'k')
end





% Timeseries of total space average
nb=W.*area;
nb(repmat(~gext.latlonmask,1,1,g.nt))=0;
Ground_delta = reshape(nansum(nansum(nb,1),2),1,[]);
nb=W.*area;
nb(repmat(gext.latlonmask,1,1,g.nt))=0;
Out_delta = reshape(nansum(nansum(nb,1),2),1,[]);
nb=W_arrival.*area;
nb(repmat(~gext.latlonmask,1,1,g.nt))=0;
Ground_arrival = reshape(nansum(nansum(nb,1),2),1,[]);
nb=W_arrival.*area;
nb(repmat(gext.latlonmask,1,1,g.nt))=0;
Out_arrival = reshape(nansum(nansum(nb,1),2),1,[]);
nb=W_departure.*area;
nb(repmat(~gext.latlonmask,1,1,g.nt))=0;
Ground_departure = reshape(nansum(nansum(nb,1),2),1,[]);
nb=W_departure.*area;
nb(repmat(gext.latlonmask,1,1,g.nt))=0;
Out_departure = reshape(nansum(nansum(nb,1),2),1,[]);

figure('position',[0 0 1000 600]);
subplot(2,1,1); hold on;legend
ylabel('Flux from the air [# bird/15min]')
plot(g.time,Ground_delta,'k','DisplayName','\Delta Departure/Arrival');
plot(g.time,Ground_arrival,'r','DisplayName','Arrival');
plot(g.time,Ground_departure,'g','DisplayName','Departure');
subplot(2,1,2); hold on;legend
plot(g.time,Out_delta,'k','DisplayName','\Delta In/Out');
plot(g.time,Out_arrival,'r','DisplayName','In');
plot(g.time,Out_departure,'g','DisplayName','Out');
ylabel('Flux from the air [# bird/15min]')

figure('position',[0 0 1000 600]);
subplot(2,1,1); hold on;legend
ylabel('In/Out the air [bird/15min]')
plot(g.time,cumsum(Ground_delta),'k','DisplayName','Flux delta Arrival/Departure');
subplot(2,1,2); hold on;legend
plot(g.time,cumsum(Out_delta),'k','DisplayName','Flux delta In/Out');



plot(g.time,-cumsum(Ground_delta),'--k');
plot(g.time,Out_delta,'r');
plot(g.time,W_allsp,'g');
xlabel('time')
ylabel('In/Out the air [bird/might]');
axis tight;



dlmwrite('W_day.dat',W_day)
dlmwrite('W_day_arrival.dat',W_day_arrival)
dlmwrite('W_day_departue.dat',W_day_departure)
dlmwrite('grid.dat',size(W_day))
dlmwrite('grid.dat',g.lat,'-append')
dlmwrite('grid.dat',g.lon,'-append')
dlmwrite('grid.dat',g.atime,'-append')

% Figure 
t=330;
figure; 
h1=subplot(2,3,1); imagesc(c(2:end-1,2:end-1,t)); colorbar;
h2=subplot(2,3,2); imagesc(vlon(2:end-1,2:end-1,t)); colorbar;
h3=subplot(2,3,3); imagesc( vlat(2:end-1,2:end-1,t)); colorbar;
h4=subplot(2,3,4); imagesc(ct-c(2:end-1,2:end-1,t)); colorbar;
h5=subplot(2,3,5); imagesc(dcdlon); colorbar;
h6=subplot(2,3,6); imagesc(dcdlat); colorbar;

linkaxes([h1,h2,h3,h4,h5,h6], 'xy');

h=figure(2);  
filename='data/SinkSourceSim';
mask_fullday = find(~reshape(all(all(isnan(W),1),2),1,[]));
worldmap([gext.lat(1) gext.lat(end)], [gext.lon(1) gext.lon(end)]); 
Frame(numel(mask_fullday)-1) = struct('cdata',[],'colormap',[]); 
colormap(brewermap([],'Spectral'))
hcoast=plotm(coastlat, coastlon,'k'); caxis([-20 20]); col=colorbar;col.Label.String = 'Source/Sink [bird/km^2/15min]';
for i_t = 1:numel(mask_fullday)

    hsurf=surfm(gext.lat2D,gext.lon2D,W(:,:,mask_fullday(i_t)));
    title(datestr(g.time(mask_fullday(i_t)))); drawnow;
%     Frame(i_t) = getframe(h);
%     [imind,cm] = rgb2ind(frame2im(Frame(i_t)),256); 
%     if i_t == 1
%         imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append');
%     end
     delete(hsurf);
end

v=VideoWriter([filename '.avi']);
v.FrameRate = 4;
v.Quality = 75;
open(v);
writeVideo(v, Frame);
close(v);

figure;
worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]); 
hcoast=plotm(coastlat, coastlon,'k');
surfm(g.lat2D(2:end-1,2:end-1),g.lon2D(2:end-1,2:end-1),nanmean(W,3));

