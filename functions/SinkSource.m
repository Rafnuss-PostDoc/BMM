% Estimation
load('data/Density_estimationMap','g')
c = (g.dens_est);
load('data/FlightSpeed_estimationMap','g')
speed = g.dens_est;
load('data/FlightDir_estimationMap','gdir')
dir = mod(gdir.dd_est,360);

% Simulation
load('data/Density_simulationMap_reassemble')
c = real_dens;
load('data/FlightSpeed_simulationMap_reassemble')
speed = real_dens;
load('data/FlightDir_simulationMap_reassemble')
dir = mod(real_dir,360);

dir=mod(90-dir,360);
vlon=speed.*cosd(dir);
vlat=speed.*sind(dir);

dlat = lldistkm([g.lat(1) 0],[g.lat(2) 0]) * ones(g.nlat-1,g.nlon)*1000;
dlon = repmat(lldistkm([g.lat(:) g.lon(1)*ones(g.nlat,1)],[g.lat(:) g.lon(2)*ones(g.nlat,1)]), 1, g.nlon-1)*1000;
dt=15*60;

dxdlon = @(x) 1/2 * ((x(2:end-1,3:end,:)-x(2:end-1,2:end-1,:))./repmat(dlon(2:end-1,2:end),1,1,g.nt) + (x(2:end-1,2:end-1,:)-x(2:end-1,1:end-2,:))./repmat(dlon(2:end-1,1:end-1),1,1,g.nt));
dxdlat = @(x) 1/2 * ((x(3:end,2:end-1,:)-x(2:end-1,2:end-1,:))./repmat(dlat(2:end,2:end-1),1,1,g.nt) + (x(2:end-1,2:end-1,:)-x(1:end-2,2:end-1,:))./repmat(dlat(1:end-1,2:end-1),1,1,g.nt));

dcdlon = dxdlon(c);
dcdlat = dxdlat(c);
dvlondlon = dxdlon(vlon);
dvlatdlat = dxdlat(vlat);

divcv = vlon(2:end-1,2:end-1,:).*dcdlon + vlat(2:end-1,2:end-1,:).*dcdlat + c(2:end-1,2:end-1,:).*( dvlondlon + dvlatdlat );
c_tp1 = c(2:end-1,2:end-1,:) - divcv.*dt;
W=c_tp1(:,:,1:end-1)-c(2:end-1,2:end-1,2:end);

W_day=nan(g.nlat,g.nlon,g.nat);
for i_t=1:g.nat-1
    idt=g.atime(i_t)+0.5<datenum(g.time) & g.atime(i_t+1)+0.5>datenum(g.time);
    W_day(2:end-1,2:end-1,i_t)=nansum(W(:,:,idt),3); 
end
W_day(~repmat(g.latlonmask,1,1,g.nat))=nan;

% save('data/SinkSourceSim.mat','W')


% Figure per night
load coastlines;
figure('position',[0 0 1000 600]); 
colormap(brewermap([],'Spectral'))
for i_t = 1:g.nat-1
    subplot(4,6,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
    surfm(g.lat2D,g.lon2D,W_day(:,:,i_t)); caxis([-300 300])
    plotm(coastlat, coastlon,'k')
end

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
mask_fullday = find(~reshape(all(all(isnan(W),1),2),g.nt-1,1));
worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]); 
Frame(numel(mask_fullday)-1) = struct('cdata',[],'colormap',[]); 
hcoast=plotm(coastlat, coastlon,'k'); caxis([-20 20]); c=colorbar;c.Label.String = 'Source/Sink [bird/km^2/15min]';
for i_t = 1:numel(mask_fullday)

    hsurf=surfm(g.lat2D(2:end-1,2:end-1),g.lon2D(2:end-1,2:end-1),W(:,:,mask_fullday(i_t)));

    title(datestr(g.time(mask_fullday(i_t)))); drawnow;
    Frame(i_t) = getframe(h);
    [imind,cm] = rgb2ind(frame2im(Frame(i_t)),256); 
    if i_t == 1
        imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append');
    end
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

