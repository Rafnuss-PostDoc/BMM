% Estimation 
load('data/Density_estimationMap','g')
load('data/Density_simulationMap_reassemble_ll')
load('data/Flight_simulationMap_real_ll')
addpath('./functions/')
load coastlines;

%% NOT GOOD YET ! 

% compute area
dlat = lldistkm([g.lat(1) g.lon(1)],[g.lat(2) g.lon(1)]);
dlon = lldistkm([g.lat2D(:,1) g.lon2D(:,1)],[g.lat2D(:,1) g.lon2D(:,2)]);
area = repmat(dlon*dlat,1,g.nlon,g.nt);
area(~g.mask_rain)=nan;

% Define the delta distance
dlat = lldistkm([g.lat(1) 0],[g.lat(2) 0]) * ones(g.nlat,g.nlon)*1000;
dlon = repmat(lldistkm([g.lat(:) g.lon(1)*ones(g.nlat,1)],[g.lat(:) g.lon(2)*ones(g.nlat,1)]), 1, g.nlon)*1000;
dt=15*60;

% Central Difference with variable \Delta x
% dxdlon = @(x) 1/2 * ((x(2:end-1,3:end,:)-x(2:end-1,2:end-1,:))./repmat(dlon(2:end-1,2:end),1,1,g.nt) + (x(2:end-1,2:end-1,:)-x(2:end-1,1:end-2,:))./repmat(dlon(2:end-1,1:end-1),1,1,g.nt));
% dxdlat = @(x) 1/2 * ((x(3:end,2:end-1,:)-x(2:end-1,2:end-1,:))./repmat(dlat(2:end,2:end-1),1,1,g.nt) + (x(2:end-1,2:end-1,:)-x(1:end-2,2:end-1,:))./repmat(dlat(1:end-1,2:end-1),1,1,g.nt));
dxdlon = @(x) 1/2 * (( [x(:,2:end,:) zeros(size(x,1),1,size(x,3))]-x)./repmat(dlon,1,1,g.nt) + (x-[zeros(size(x,1),1,size(x,3)) x(:,1:end-1,:)])./repmat(dlon,1,1,g.nt));
dxdlat = @(x) 1/2 * (( [x(2:end,:,:) zeros(1,size(x,2),size(x,3))]-x)./repmat(dlat,1,1,g.nt) + (x-[zeros(1,size(x,2),size(x,3)) x(1:end-1,:,:)])./repmat(dlat,1,1,g.nt));



vlon = zeros(g.nlat,g.nlon, g.nt);
vlat = zeros(g.nlat,g.nlon, g.nt);
c = zeros(g.nlat,g.nlon, g.nt);

latlonmask = repmat(g.latlonmask,1,1,g.nt);


nb_real = 10;%size(real_dens_ll,3);
W = nan(g.nlat,g.nlon,g.nt,nb_real);

 
for i=1:nb_real
    vlon(latlonmask) = real_u_ll(:,:,i);
    vlat(latlonmask) = real_v_ll(:,:,i);
    c(latlonmask)    = real_dens_ll(:,:,i);
    
    dcvdlon = dxdlon(c.*vlon);
    dcvdlon2D = dcvdlon;
    dcvdlon2D(:,2:end-1,:) = 1/4*( dcvdlon(:,1:end-2,:) +2*dcvdlon(:,2:end-1,:) +dcvdlon(:,3:end,:) ); 
    dcvdlat = dxdlon(c.*vlat);
    dcvdlat2D = dcvdlat;
    dcvdlat2D(2:end-1,:,:) = 1/4*( dcvdlat(1:end-2,:,:) +2*dcvdlat(2:end-1,:,:) +dcvdlat(3:end,:,:) ); 
    F = dcvdlon2D+dcvdlat2D;
    W(:,:,1:end-1,i) = c(:,:,2:end) - c(:,:,1:end-1) + 0.5*dt*(F(:,:,2:end)+F(:,:,1:end-1));
    i
end
% c per unit volume







% Figure for all time
W_all=nan(g.nlat,g.nlon,nb_real);
for i=1:nb_real
    W_all(:,:,i)=nansum(W(:,:,:,i),3);
end
W_all(repmat(~g.latlonmask,1,1,nb_real))=nan;

figure('position',[0 0 1000 600]); 
colormap(brewermap([],'Spectral'))
for i=1:nb_real
    subplot(2,(nb_real)/2,i)
    worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
    surfm(g.lat2D,g.lon2D,W_all(:,:,i)); 
    caxis([min(W_all(:)) max(W_all(:))])
    col=colorbar;col.Label.String = 'Bird Sink/source [bird/km^2]';
    plotm(coastlat, coastlon,'k')
end
subplot(2,(nb_real+1)/2,nb_real+1)
worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
surfm(g.lat2D,g.lon2D,mean(W_all,3)); 
caxis([min(W_all(:)) max(W_all(:))])
col=colorbar;col.Label.String = 'Bird Sink/source [bird/km^2]';
plotm(coastlat, coastlon,'k')

% Figure for all space
W_alls=nan(g.nt,nb_real);
for i=1:nb_real
    W_alls(:,i)=nansum(nansum(W(:,:,:,i).*area,1),2);
end
figure('position',[0 0 1000 600]); 
plot(g.time,W_alls)
xlabel('time')
ylabel('In/Out the air [bird]')


W_day=nan(g.nlat,g.nlon,g.nat,nb_real);
for i_t=1:g.nat-1
    idt=g.atime(i_t)-0.5<datenum(g.time) & g.atime(i_t)+0.5>datenum(g.time);
    for i=1:nb_real
        W_day(:,:,i_t,i)=nansum(W(:,:,idt,i),3);
    end
end
W_day(~repmat(g.latlonmask,1,1,g.nat,nb_real))=nan;

% save('data/SinkSourceSim.mat','W')


% Figure per night
load coastlines;
for i=1:nb_real
    figure('position',[0 0 1000 600]); 
    %colormap(brewermap([],'Spectral'))
    for i_t = 1:g.nat-1
        subplot(4,6,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
        surfm(g.lat2D,g.lon2D,W_day(:,:,i_t,i)); caxis([-100 100])
        plotm(coastlat, coastlon,'k')
    end
end
figure('position',[0 0 1000 600]); 
colormap(brewermap([],'Spectral'))
for i_t = 1:g.nat-1
    subplot(4,6,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
    surfm(g.lat2D,g.lon2D,mean(W_day(:,:,i_t,:),4)); caxis([-300 300])
    plotm(coastlat, coastlon,'k')
end


dlmwrite('W_day.dat',W_day)
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

