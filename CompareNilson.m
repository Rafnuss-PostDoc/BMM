
%% At radar location
% Loading
load('./data/dc_corr.mat');
VP=readtable('./data/vp_processed_70_radars_20160919_20161009.csv');


% Compare data to data
for i_r=1:numel(dc)
    id = strcmp(dc(i_r).name,VP.radar_id);
    if sum(id)>0
        % duplicate for 'searl' 2016-10-09 00:00:00 
        try
            dc(i_r).dens2=reshape(VP.dens(id),size(dc(i_r).dens'))';
            dc(i_r).ff2=reshape(VP.ff(id),size(dc(i_r).ff'))';
            dc(i_r).dd2=reshape(VP.dd(id),size(dc(i_r).dd'))';
        catch
            [C,ia,ic] =unique(VP(id,[2 3]));
            id=find(id);
            dc(i_r).dens2=reshape(VP.dens(id(ia)),size(dc(i_r).dens'))';
            dc(i_r).ff2=reshape(VP.ff(id(ia)),size(dc(i_r).ff'))';
            dc(i_r).dd2=reshape(VP.dd(id(ia)),size(dc(i_r).dd'))';
        end
    end
end

% VP export file
for i_d=1:numel(dc)
    if dc(i_d).interval>0
        figure(6); clf; set(gcf, 'Color', 'w');

        subplot(2,1,1); hold on;
        imagesc(datenum(dc(i_d).time), dc(i_d).interval*(1/2:double(dc(i_d).levels)), (dc(i_d).ff)','AlphaData',~isnan(dc(i_d).dens'))
        xlabel('date'); ylabel('[m]'); title('Density me [log10 bird/m^3]');colorbar; caxis([0 30])
        axis([datenum(start_date) datenum(end_date) 0 5000]); datetick('x','dd/mm','keeplimits'); set(gca, 'YDir', 'normal')
        
        subplot(2,1,2); hold on;
        imagesc(datenum(dc(i_d).time), dc(i_d).interval*(1/2:double(dc(i_d).levels)), (dc(i_d).ff2)','AlphaData',~isnan(dc(i_d).dens'))
        plot([datenum(start_date) datenum(end_date)],[dc(i_d).height dc(i_d).height],'r','linewidth',2); caxis([0 30])
        ylabel('[m]'); title('Density Nilson [log10 bird/m^3]');datetick('x','dd/mm','keeplimits');   colorbar;
        datetick('x'); axis([datenum(start_date) datenum(end_date) 0 5000]); set(gca, 'YDir', 'normal');
        export_fig(['figure\comparisonNilson\' num2str(i_d) '_' dc(i_d).name '_ff.png'],'-png')
    end
end




% Nilson's way
uniqueDate = round(datenum(datetime(start_date-1:end_date-1)));
dd_avg=nan(numel(dc), numel(uniqueDate));
for i_r=1:numel(dc)
    %deltalat=dc(i_r).ff * 60*60/1000 .* cosd(dc(i_r).dd);
    %deltalon=dc(i_r).ff * 60*60/1000 .* sind(dc(i_r).dd);
    try
        dc(i_r).scoret = datenum(dc(i_r).time'-mean([dc(i_r).sunrise;dc(i_r).sunset])) ./ datenum(dc(i_r).sunrise-dc(i_r).sunset)*2;
        idx = dc(i_r).scoret(:)>=-1 & dc(i_r).scoret(:)<=1 & dc(i_r).dens2>5;

        for i_d=1:numel(uniqueDate)
            idx2=round(datenum(dc(i_r).sunset))==uniqueDate(i_d);
            idx3=idx;
            idx3(~idx2,:)=0;
            a = dc(i_r).dd2(idx3);
            a = a(~isnan(a));
            %             vx = sind(dc(i_r).dd2(idx3));
            %             vy = cosd(dc(i_r).dd2(idx3));
            %             dd_avg(i_r,i_d) = atand(nanmean(vx)/nanmean(vy))+180;
            if ~isempty(a)
                dd_avg(i_r,i_d) = rad2deg(circ_mean(deg2rad(a)-pi)+pi);
            end
            % dd_avg(i_r,i_d) = nansum(dc(i_r).dd(idx3) .* dc(i_r).dens(idx3)) / nansum(dc(i_r).dens(idx3));
        end
    end
end

figure('position',[0 0 1000 600]); hold on;
worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
plotm(coastlat, coastlon,'k')
s=.1;
for i_r=1:numel(dc)
    quiverm(dc(i_r).lat,dc(i_r).lon,cosd(nanmean(dd_avg(i_r,:))),sind(nanmean(dd_avg(i_r,:))))
end



% My way

for i_r=1:numel(dc)
    idx = data.i_r==i_r;
    alpha = deg2rad(data.dd(idx));
    % w_MTR = data.ff(idx).* data.dens(idx) ./ nansum(data.ff(idx).* data.dens(idx));
    w_MTR = data.dens(idx) ./ nansum(data.dens(idx));
    r = nansum(w_MTR.*exp(1i*alpha));
    dd_avg(i_r) = rad2deg(atan2(imag(r), real(r)));
end

figure('position',[0 0 1000 600]); hold on;
worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
plotm(coastlat, coastlon,'k')
s=.1;
for i_r=1:numel(dc)
    quiverm(dc(i_r).lat,dc(i_r).lon,cosd(dd_avg(i_r)),sind(dd_avg(i_r)))
end


%% Full map

% Load data from BMM
load('data/Density_estimationMap','g')
c = (g.dens_est); % bird/km^2

load('data/FlightSpeed_estimationMap','g')
speed = g.dens_est * 60*60/1000;% convert to km/hr

load('data/FlightDir_estimationMap','gdir')
dir = mod(gdir.dd_est,360); % deg

load coastlines;

% MTR
MTR = c .* speed;

% Flight Quiver
dir=mod(90-dir,360);
deltalat=speed.*cosd(dir);
deltalon=speed.*sind(dir);





% MTR over all period
figure('position',[0 0 1000 600]); hold on;
worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
surfm(g.lat2D,g.lon2D,nanmean(MTR,3))
contourm(g.lat2D,g.lon2D,nanmean(MTR,3),0:250:1250,'k');
plotm(coastlat, coastlon,'k')
circlem([dc.lat],[dc.lon],10,'filled')
deltalat_avg = nansum(deltalat.*c,3)./nansum(c,3);
deltalon_avg = nansum(deltalon.*c,3)./nansum(c,3);
d=5;
quiverm(g.lat2D(1:d:end,1:d:end),g.lon2D(1:d:end,1:d:end),deltalat_avg(1:d:end,1:d:end),deltalon_avg(1:d:end,1:d:end),'k')


% MTR per night
MTR_day=nan(g.nlat,g.nlon,g.nat);
for i_t=1:g.nat-1
    idt=g.atime(i_t)-0.5<datenum(g.time) & g.atime(i_t)+0.5>datenum(g.time);
    MTR_day(:,:,i_t)=nanmean(MTR(:,:,idt),3); 
end
MTR_day(~repmat(g.latlonmask,1,1,g.nat))=nan;


figure('position',[0 0 1000 600]); 
for i_t = 1:g.nat-1
    subplot(4,6,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
    surfm(g.lat2D,g.lon2D,MTR_day(:,:,i_t)); caxis([0 3000])
    plotm(coastlat, coastlon,'k')
end




% Average of the period of time
dates=datetime({'16-Sep-2016 12:00','23-Sep-2016 12:00','30-Sep-2016 12:00','7-Oct-2016 12:00','23-Oct-2016 12:00'},'InputFormat','dd-MMM-yyy HH:mm');
MTR_week=nan(g.nlat,g.nlon,numel(dates)-1);
for i_t=1:numel(dates)-1
    idt=datenum(dates(i_t))<g.atime & datenum(dates(i_t+1))>g.atime;
    MTR_week(:,:,i_t)=nanmean(MTR_day(:,:,idt),3); 
end

load coastlines;
figure('position',[0 0 1000 600]); 
for i_t = 1:numel(dates)-1
    subplot(1,4,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
    surfm(g.lat2D,g.lon2D,MTR_week(:,:,i_t)); caxis([0 3000])
    plotm(coastlat, coastlon,'k')
end