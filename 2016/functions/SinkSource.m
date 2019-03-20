% Estimation 
load('data/Density_estimationMap','gd','g')
load('data/Flight_estimationMap','guv')
addpath('./functions/')
load coastlines;

vlon= guv.u_est; % m/s (+) east, (-) west
vlat=guv.v_est; % m/s (+) north, (-) south
c = gd.dens_est;





% compute area
dlat = lldistkm([g.lat(1) g.lon(1)],[g.lat(2) g.lon(1)]);
dlon = lldistkm([g.lat2D(:,1) g.lon2D(:,1)],[g.lat2D(:,1) g.lon2D(:,2)]);
area = repmat(dlon*dlat,1,g.nlon,g.nt);
area(~g.mask_rain)=nan;

% Define the delta distance
dlat = lldistkm([g.lat(1) 0],[g.lat(2) 0]) * ones(g.nlat-1,g.nlon)*1000;
dlon = repmat(lldistkm([g.lat(:) g.lon(1)*ones(g.nlat,1)],[g.lat(:) g.lon(2)*ones(g.nlat,1)]), 1, g.nlon-1)*1000;
dt=15*60;

% Central Difference with variable \Delta x
dxdlon = @(x) 1/2 * ((x(2:end-1,3:end,:)-x(2:end-1,2:end-1,:))./repmat(dlon(2:end-1,2:end),1,1,g.nt) + (x(2:end-1,2:end-1,:)-x(2:end-1,1:end-2,:))./repmat(dlon(2:end-1,1:end-1),1,1,g.nt));
dxdlat = @(x) 1/2 * ((x(3:end,2:end-1,:)-x(2:end-1,2:end-1,:))./repmat(dlat(2:end,2:end-1),1,1,g.nt) + (x(2:end-1,2:end-1,:)-x(1:end-2,2:end-1,:))./repmat(dlat(1:end-1,2:end-1),1,1,g.nt));



% Option 1
% dcdlon = dxdlon(c);
% dcdlat = dxdlat(c);
% dvlondlon = dxdlon(vlon);
% dvlatdlat = dxdlat(vlat);
% 
% F = vlon(2:end-1,2:end-1,:).*dcdlon + vlat(2:end-1,2:end-1,:).*dcdlat + c(2:end-1,2:end-1,:).*( dvlondlon + dvlatdlat );
% c_tp1 = c(2:end-1,2:end-1,:) - F.*dt;
% W2=c_tp1(:,:,1:end-1)-c(2:end-1,2:end-1,2:end);


W = nan(g.nlat,g.nlon,g.nt);

dcvdlon = dxdlon(c.*vlon);
dcvdlon2D = dcvdlon;
dcvdlon2D(:,2:end-1,:) = 1/4*( dcvdlon(:,1:end-2,:) +2*dcvdlon(:,2:end-1,:) +dcvdlon(:,3:end,:) ); 
dcvdlat = dxdlon(c.*vlat);
dcvdlat2D = dcvdlat;
dcvdlat2D(2:end-1,:,:) = 1/4*( dcvdlat(1:end-2,:,:) +2*dcvdlat(2:end-1,:,:) +dcvdlat(3:end,:,:) ); 
F = dcvdlon+dcvdlat;
c_tp1 = c(2:end-1,2:end-1,1:end-1) - 0.5*dt*(F(:,:,2:end)+F(:,:,1:end-1));
W(2:end-1,2:end-1,1:end-1) = c_tp1-c(2:end-1,2:end-1,2:end);
W=-W;

% save('data/SinkSourceSim.mat','W')






% Compute flux out
fluxlon = c .* vlon .* repmat(dlon(:,1),1,g.nlon,g.nt) * dt;
fluxlat = c .* vlat .* repmat(dlat(1),g.nlat,g.lon,g.nt) * dt;

% Lat flux
fluxlonfm = fillmissing(fluxlon,'linear');
Fg = griddedInterpolant(g.lat3D,g.lon3D,datenum(g.time3D),fluxlon);
fluxlatdlat = Fg( repmat(g.lat(1:end-1)'+diff(g.lat')/2,1,g.nlon,g.nt) , g.lon3D(1:end-1,:,:) , datenum(g.time3D(1:end-1,:,:)) );

latborder = diff(g.latlonmask,1,1);

fluxlatdlat(repmat(latborder,1,1,g.nt)>0)

for i=1:g.nt
    imagesc(fluxlatdlat(:,:,i));
    drawnow
end

latfluxlat = repmat(g.lat(1:end-1)'+diff(g.lat')/2,1,g.nlon);
latfluxlon = g.lat2D(1:end-1,:);
latflux = Fg(repmat(latfluxlat(~latborder==0),1,g.nt),repmat(latfluxlon(~latborder==0),1,g.nt),repmat(datenum(g.time),sum(~latborder(:)==0),1));


lonborder = diff(g.latlonmask,1,2);

fluxlonout = repmat(lonborder,1,1,g.nt) .* fluxlon;
fluxlatout = repmat(latborder,1,1,g.nt) .* fluxlat;



















% Positiv and negativ flux
Wm=W; Wp=W;
Wm(Wm>0)=0;
Wp(Wp<0)=0;





% Map of total time average
W_all=nansum(W(:,:,:),3);
W_all(repmat(~g.latlonmask,1,1))=nan;

figure('position',[0 0 1000 600]); 
colormap(brewermap([],'Spectral'))
worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
surfm(g.lat2D,g.lon2D,W_all); 
caxis([min(W_all(:)) max(W_all(:))])
c=colorbar;c.Label.String = 'Bird Sink/source [bird/km^2]';
plotm(coastlat, coastlon,'k')



% Maps of daily average
W_day=nan(g.nlat,g.nlon,g.nat);
W_dayp=nan(g.nlat,g.nlon,g.nat);
W_daym=nan(g.nlat,g.nlon,g.nat);
for i_t=1:g.nat-1
    idt=g.atime(i_t)-0.5<datenum(g.time) & g.atime(i_t)+0.5>datenum(g.time);
    W_day(:,:,i_t)=nansum(W(:,:,idt),3);
    W_dayp(:,:,i_t)=nansum(Wp(:,:,idt),3);
    W_daym(:,:,i_t)=nansum(Wm(:,:,idt),3);
end
W_day(~repmat(g.latlonmask,1,1,g.nat,nb_real))=nan;
W_dayp(~repmat(g.latlonmask,1,1,g.nat,nb_real))=nan;
W_daym(~repmat(g.latlonmask,1,1,g.nat,nb_real))=nan;

figure('position',[0 0 1000 600]); 
colormap(brewermap([],'Spectral'))
for i_t = 1:g.nat-1
    subplot(4,6,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
    surfm(g.lat2D,g.lon2D,W_day(:,:,i_t)); caxis([-100 100])
    plotm(coastlat, coastlon,'k')
end

figure('position',[0 0 1000 600]); 
for i_t = 1:g.nat-1
    subplot(4,6,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
    surfm(g.lat2D,g.lon2D,W_dayp(:,:,i_t)); caxis([0 100])
    plotm(coastlat, coastlon,'k')
end

figure('position',[0 0 1000 600]); 
for i_t = 1:g.nat-1
    subplot(4,6,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
    surfm(g.lat2D,g.lon2D,W_daym(:,:,i_t)); caxis([-100 0])
    plotm(coastlat, coastlon,'k')
end


% Timeseries of total space average
W_alls=reshape(nansum(nansum(W.*area,1),2),1,[]);
W_allsm=reshape(nansum(nansum(Wm.*area,1),2),1,[]);
W_allsp=reshape(nansum(nansum(Wp.*area,1),2),1,[]);

figure('position',[0 0 1000 600]); hold on;
ylabel('In/Out the air [bird/15min]')
plot(g.time,W_alls,'k');
plot(g.time,-cumsum(W_alls),'--k');
plot(g.time,W_allsm,'r');
plot(g.time,W_allsp,'g');
xlabel('time')
ylabel('In/Out the air [bird/might]');
axis tight;



dlmwrite('W_day.dat',W_day)
dlmwrite('W_day_arrival.dat',W_daym)
dlmwrite('W_day_departue.dat',W_dayp)
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
worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]); 
Frame(numel(mask_fullday)-1) = struct('cdata',[],'colormap',[]); 
colormap(brewermap([],'Spectral'))
hcoast=plotm(coastlat, coastlon,'k'); caxis([-10 10]); col=colorbar;col.Label.String = 'Source/Sink [bird/km^2/15min]';
for i_t = 1:numel(mask_fullday)

    hsurf=surfm(g.lat2D(2:end-1,2:end-1),g.lon2D(2:end-1,2:end-1),W(:,:,mask_fullday(i_t)));
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

