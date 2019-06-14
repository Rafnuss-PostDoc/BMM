%% Sink/source version 

% Estimation 
load('data/Density_estimationMap','gd','g')
load('data/Flight_estimationMap','guv')
addpath('./functions/')
load coastlines;

dn = 1;
gext.lat = fillmissing([nan(1,dn) g.lat nan(1,dn)],'linear');
gext.lon = fillmissing([nan(1,dn) g.lon nan(1,dn)],'linear');
gext.nlat = numel(gext.lat);
gext.nlon = numel(gext.lon);
[gext.lon2D, gext.lat2D] = meshgrid(gext.lon,gext.lat);

% compute area
dlat2 = lldistkm([g.lat(1) g.lon(1)],[g.lat(2) g.lon(1)]);
dlon2 = lldistkm([g.lat2D(:,1) g.lon2D(:,1)],[g.lat2D(:,1) g.lon2D(:,2)]);
area = repmat(dlon2*dlat2,1,g.nlon,g.nt-1);
dlatext2 = lldistkm([gext.lat(1) gext.lon(1)],[gext.lat(2) gext.lon(1)]);
dlonext2 = lldistkm([gext.lat2D(:,1) gext.lon2D(:,1)],[gext.lat2D(:,1) gext.lon2D(:,2)]);
areaext = repmat(dlonext2*dlatext2,1,gext.nlon,g.nt-1);


% Define the delta distance
dlat = lldistkm([g.lat(1) 0],[g.lat(2) 0]) * ones(g.nlat,g.nlon); % ° -> km
dlon = repmat(lldistkm([g.lat(:) g.lon(1) * ones(g.nlat,1)],[g.lat(:) g.lon(2)*ones(g.nlat,1)]), 1, g.nlon); % ° -> km
dlatext = lldistkm([gext.lat(1) 0],[gext.lat(2) 0]) * ones(gext.nlat,gext.nlon); % ° -> km
dlonext = repmat(lldistkm([gext.lat(:) gext.lon(1) * ones(gext.nlat,1)],[gext.lat(:) gext.lon(2)*ones(gext.nlat,1)]), 1, gext.nlon); % ° -> km


dt=15/60; % -> 15 min -> 15/60 hr

vlat = guv.v_est/1000*60*60; % m/s -> km/h (+) north, (-) south
vlon = guv.u_est/1000*60*60; % m/s -> km/h (+) east, (-) west
rho =  gd.dens_est; % bird/km^2
rho(repmat(g.latlonmask,1,1,g.nt) & isnan(rho))=0;

% add a nan layer in lat or lon direction respectively for the computation
% of the flux at 1/2
Philat_pad = padarray(rho.*vlat,[1 0 0],nan);
Philon_pad = padarray(rho.*vlon,[0 1 0],nan);

% Compute the i-1/2 and i+1/2 flux .
% compute even if either previous or next cells is not nan. (boundary
% fluxes)
Philat_h = movmean(Philat_pad,[0 1],1,'omitnan','Endpoints','discard'); 
Philon_h = movmean(Philon_pad,[0 1],2,'omitnan','Endpoints','discard');

% Compute the delta flux / delta distance
% diff elimiate if either i-1/2 or i+1/2 flux is nan.
% -> not usefull for flux on boundary
dPhilatdlat = diff(Philat_h,1,1)./ repmat(dlat,1,1,g.nt);
dPhilondlon = diff(Philon_h,1,2)./repmat(dlon,1,1,g.nt);

% perform a average over the perpendicular axes to account for ... ?
% dcvdlonP = padarray(dPhilondlon,[1 1],nan);
% dcvdlatP = padarray(dPhilatdlat,[1 1],nan);
% dcvdlon2D = nanmean(cat(4,dcvdlonP(2:end-1,1:end-2,:), dcvdlonP(2:end-1,2:end-1,:), dcvdlonP(2:end-1,2:end-1,:), dcvdlonP(2:end-1,3:end,:)),4);
% dcvdlat2D = nanmean(cat(4,dcvdlatP(1:end-2,2:end-1,:), dcvdlatP(1:end-2,2:end-1,:), dcvdlatP(1:end-2,2:end-1,:), dcvdlatP(3:end,2:end-1,:)),4);

F = dPhilatdlat + dPhilondlon; 
F_h = movmean(F,[0 1],3,'omitnan','Endpoints','discard');
Win = diff(rho,1,3) + dt*F_h;
% bird/km^2 (*15min/15min)


% Check for a step
i=20; i=70;
figure;
subplot(3,2,1); imagesc(Philon_pad(:,:,i)); set(gca,'ydir','normal'); colorbar;
subplot(3,2,2); imagesc(Philat_pad(:,:,i)); set(gca,'ydir','normal'); colorbar;
subplot(3,2,3); imagesc(dPhilondlon(:,:,i)); set(gca,'ydir','normal'); colorbar;
subplot(3,2,4); imagesc(dPhilatdlat(:,:,i)); set(gca,'ydir','normal'); colorbar;
subplot(3,2,5); imagesc(rho(:,:,i+1) - rho(:,:,i)); set(gca,'ydir','normal'); colorbar;
subplot(3,2,6); imagesc(F_h(:,:,i)); set(gca,'ydir','normal'); colorbar;



% Inside / Outside fluxes
% add 0 padding for the outer zone
Philat_h_0=padarray(Philat_h,[1 1 0],0);
% replace nan by zero for the diff to work
Philat_h_0(isnan(Philat_h_0))=0;
% Flux [bird/km^2*km/hr]=[bird/km/hr] -> NB [bird/hr] by multiplying by
% cell size perepndicular (lat-> lon)
Nlat = diff(Philat_h_0,1,1).* repmat(dlonext,1,1,g.nt);
Philon_h_0=padarray(Philon_h,[1 1 0],0);
Philon_h_0(isnan(Philon_h_0))=0;
Nlon = diff(Philon_h_0,1,2)./ repmat(dlatext,1,1,g.nt);

% [bird/hr] -> [bird]
NBout = dt*(Nlat+Nlon); 
% remove value inside the domain
NBout( repmat(padarray(g.latlonmask,[1 1 0],0),1,1,g.nt) )=0;
% average for number of bird between two time step.
NBout = movmean(NBout,[0 1],3,'omitnan','Endpoints','discard');



% Total W
Wint=padarray(Win,[1 1 0],0);
Wint(isnan(Wint))=0;
Woutt = NBout ./ areaext;
W=Wint+Woutt;
W(W==0)=nan;


% save('data/SinkSourceSim.mat','W')




% Map of total time average
W_all=nansum(Win(:,:,:),3);
W_all(W_all==0)=NaN;


figure('position',[0 0 1000 600]); 
colormap(brewermap([],'Spectral'))
worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow('worldrivers.shp','Color', 'blue')
surfm(g.lat2D,g.lon2D,W_all); 
caxis([min(W_all(:)) max(W_all(:))])
col=colorbar;col.Label.String = 'Bird Sink/source [bird/km^2]';
plotm(coastlat, coastlon,'k')



% Positiv and negativ flux
W_arrival=Win; W_departure=Win;
W_arrival(W_arrival>=0)=0;
W_departure(W_departure<=0)=0;


% Maps of daily average
W_day=nan(g.nlat,g.nlon,g.nat);
W_day_departure=nan(g.nlat,g.nlon,g.nat);
W_day_arrival=nan(g.nlat,g.nlon,g.nat);
for i_t=1:g.nat-1
    idt=g.atime(i_t)-0.5<datenum(g.time) & g.atime(i_t)+0.5>datenum(g.time);
    W_day(:,:,i_t)=nansum(Win(:,:,idt),3);
    W_day_departure(:,:,i_t)=nansum(W_departure(:,:,idt),3);
    W_day_arrival(:,:,i_t)=nansum(W_arrival(:,:,idt),3);
end
W_day(W_day==0)=nan;
W_day_departure(W_day_departure==0)=nan;
W_day_arrival(W_day_arrival==0)=nan;

figure('position',[0 0 1000 600]); 
%colormap(brewermap([],'Spectral'))
for i_t = 1:g.nat-1
    subplot(4,6,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
    surfm(g.lat2D,g.lon2D,W_day(:,:,i_t)); caxis([-100 100])
    plotm(coastlat, coastlon,'k')
end

figure('position',[0 0 1000 600]); 
for i_t = 1:g.nat-1
    subplot(4,6,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
    surfm(g.lat2D,g.lon2D,W_day_departure(:,:,i_t)); caxis([0 200])
    plotm(coastlat, coastlon,'k')
end

figure('position',[0 0 1000 600]); 
for i_t = 1:g.nat-1
    subplot(4,6,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
    surfm(g.lat2D,g.lon2D,W_day_arrival(:,:,i_t)); colorbar; caxis([-200 0])
    plotm(coastlat, coastlon,'k')
end




% Positiv and negativ flux
NBout_out=NBout; NBout_in=NBout;
NBout_out(NBout_out>=0)=0;
NBout_in(NBout_in<=0)=0;


% Timeseries of total space average
nb=Win.*area;
In_delta = reshape(nansum(nansum(nb,1),2),1,[]);
nb=W_arrival.*area;
In_arrival = reshape(nansum(nansum(nb,1),2),1,[]);
nb=W_departure.*area;
In_departure = reshape(nansum(nansum(nb,1),2),1,[]);

Out_delta=reshape(nansum(nansum(NBout,1),2),1,[]);
Out_out=reshape(nansum(nansum(NBout_out,1),2),1,[]);
Out_in=reshape(nansum(nansum(NBout_in,1),2),1,[]);



figure('position',[0 0 1000 600]);
subplot(2,1,1); hold on;legend
ylabel('Flux from the air [# bird]')
plot(g.time(1:end-1),In_delta,'k','DisplayName','\Delta Departure/Arrival');
plot(g.time(1:end-1),In_arrival,'r','DisplayName','Arrival');
plot(g.time(1:end-1),In_departure,'g','DisplayName','Departure');
subplot(2,1,2); hold on;legend
plot(g.time(1:end-1),Out_delta,'k','DisplayName','\Delta In/Out');
plot(g.time(1:end-1),Out_in,'r','DisplayName','In');
plot(g.time(1:end-1),Out_out,'g','DisplayName','Out');
ylabel('Flux from the air [# bird]')

figure('position',[0 0 1000 600]);
subplot(2,1,1); hold on;legend
ylabel('In/Out the air [bird]')
plot(g.time(1:end-1),cumsum(In_delta),'k','DisplayName','Flux delta Arrival/Departure');
subplot(2,1,2); hold on;legend
plot(g.time(1:end-1),cumsum(Out_delta),'k','DisplayName','Flux delta In/Out');



plot(g.time,-cumsum(In_delta),'--k');
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
h1=subplot(2,3,1); imagesc(rho(2:end-1,2:end-1,t)); colorbar;
h2=subplot(2,3,2); imagesc(vlon(2:end-1,2:end-1,t)); colorbar;
h3=subplot(2,3,3); imagesc( vlat(2:end-1,2:end-1,t)); colorbar;
h4=subplot(2,3,4); imagesc(ct-rho(2:end-1,2:end-1,t)); colorbar;
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
keyboard
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
surfm(g.lat2D(2:end-1,2:end-1),g.lon2D(2:end-1,2:end-1),nanmean(Win,3));

