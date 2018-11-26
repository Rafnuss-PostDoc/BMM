% Load data frp, EBP
R = readtable('data\EBP_files\Records_2016.csv');
E = readtable('data\EBP_files\Events_2016.csv');

% Load data from BMM
load('data/Density_estimationMap','g')
c = (g.dens_est); % bird/km^2

load('data/FlightSpeed_estimationMap','g')
speed = g.dens_est * 60*60/1000; % convert to km/hr


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
    subplot(1,4,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
    surfm(g.lat2D,g.lon2D,MTR_week(:,:,i_t)); caxis([0 3000])
    plotm(coastlat, coastlon,'k')
end




% Record and Event restructure
RE.x = sort(unique(R.MAPX));
RE.y = sort(unique(R.MAPY));
RE.t = sort(unique(R.DATE_WEEK));
[RE.X, RE.Y, RE.T] = meshgrid(RE.x,RE.y,RE.t);

vn = E.Properties.VariableNames(7:end);
for i_vn=1:numel(vn)
    RE.(string(vn(i_vn)))= nan(size(RE.X));
end

for e=1:size(E,1)
    id = RE.X==E.MAPX(e) & RE.Y==E.MAPY(e) & RE.T==E.DATE_WEEK(e);
    for i_vn=1:numel(vn)
        RE.(string(vn(i_vn)))(id) = E.(char(vn(i_vn)))(e);
    end
end

RE.species = sort(unique(R.SPECIES_CD));
vn = R.Properties.VariableNames(8:end);
for i_sp=1:numel(RE.species)
    RE.(string(RE.species(i_sp)))=struct();
    for i_vn=1:numel(vn)
        RE.(string(RE.species(i_sp))).(string(vn(i_vn))) = nan(size(RE.X));
    end
end

for r=1:size(R,1)
    id = RE.X==R.MAPX(r) & RE.Y==R.MAPY(r) & RE.T==R.DATE_WEEK(r);
    for i_sp=1:numel(RE.species)
        for i_vn=1:numel(vn)
            RE.(string(RE.species(i_sp))).(string(vn(i_vn)))(id) = R.(char(vn(i_vn)))(r);
        end
    end
end

% Convert coordinate
[lat,lon] = projinv(proj,x,y)

['http://geodesy.geo.admin.ch/reframe/lv95towgs84?easting=' num2str(RE.x(end/2)) '&northing=' num2str(RE.y(end/2))]


% Test CRO

RE.OENOEN.RECORDS



for i=1:5
    subplot(1,5,i); imagesc()
    subplot(1,4,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
    surfm(RE.X,RE.Y,RE.OENOEN.RECORDS(:,:,i));
    plotm(coastlat, coastlon,'k')
end

























