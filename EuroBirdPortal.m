% Load data frp, EBP
R = readtable('data\EBP_files\Records_2016.csv');
E = readtable('data\EBP_files\Events_2016.csv');

% Load data from BMM
load('data/Density_estimationMap','g')
c = (g.dens_est);

load('data/FlightSpeed_estimationMap','g')
speed = g.dens_est * 60*60/1000;

% load('data/SinkSource.mat')




% MTR per night
MTR = c .* speed;
MTR_day=nan(g.nlat,g.nlon,g.nat);
for i_t=1:g.nat-1
    idt=g.atime(i_t)-0.5<datenum(g.time) & g.atime(i_t)+0.5>datenum(g.time);
    MTR_day(:,:,i_t)=nanmean(MTR(:,:,idt),3); 
end
MTR_day(~repmat(g.latlonmask,1,1,g.nat))=nan;

load coastlines;
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
    subplot(2,2,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
    surfm(g.lat2D,g.lon2D,MTR_week(:,:,i_t)); caxis([0 3000])
    plotm(coastlat, coastlon,'k')
end
