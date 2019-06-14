addpath('functions/')

%% 1. Download value
start_date='01-Jan-2018 00:00:00';
end_date='01-Jan-2019 00:00:00';
quantity = {'dens','ff','dd','DBZH','eta','sd_vvp'};

d = download_vp(start_date,end_date, quantity);

save('data/d2018all.mat','d','start_date','end_date','quantity','-v7.3');

% Figure
for i_d=1:numel(d)
    clf;
    subplot(2,1,1); hold on;
    surf(datenum(d(i_d).time), d(i_d).interval/1000*(1/2:double(d(i_d).levels)), (d(i_d).eta)','EdgeColor','none')
    xlabel('Date'); ylabel('Altitude [km]'); c=colorbar; c.Label.String='Raw reflectivity'; %caxis([0 5])
    datetick('x','dd mmm'); axis([datenum(start_date) datenum(end_date)-1 0 5]); view(2)
    
    subplot(2,1,2); hold on;
    surf(datenum(d(i_d).time), d(i_d).interval/1000*(1/2:double(d(i_d).levels)), log10(d(i_d).dens+1)','EdgeColor','none')
    rectangle('Position',[0 -1 datenum(end_date)+2 1+d(i_d).height/1000], 'FaceColor',[1 1 1],'EdgeColor','k','LineWidth',2)
    xlabel('Date'); ylabel('Altitude [km]'); c=colorbar; c.Label.String='Bird Density [Log bird/km^3]'; caxis([0 2]);%max(log10(d(i_d).dens(:)+1))])
    datetick('x','dd mmm'); axis([datenum(start_date) datenum(end_date)-1 0 5]); set(gca, 'YDir', 'normal')

    drawnow
    print(['figure/' num2str(i_d) ,'_', d(i_d).name '.png'],'-dpng', '-r300')
end


%% 2. Delete entire radar 
load('data/d2018all.mat'); addpath('functions/')

% Observed data
% CleaningV(d)

% Remove radars with wrong, imcomplete data: 76 radars overs 123
C=readtable('data/Cleaning.xlsx');


for i_d=numel(d):-1:1
    id = strcmp(C.Name,d(i_d).name);
    if ~C.Keep(id)
        d(i_d)=[];
    else
        d(i_d).scatter_lim=C.FirstBinOkInEta(id);
    end
end


% Compute number of data left
nb_date=zeros(numel(d),365);
for i_d=1:numel(d)
    tmp = tabulate(round(datenum(d(i_d).time))-datenum('1-1-2018'));
    nb_date(i_d,tmp(:,1))= tmp(:,2);
end

% figure; hold on;
% worldmap([min([d.lat]) max([d.lat])], [min([d.lon]) max([d.lon])]);  
% geoshow( 'landareas.shp', 'FaceColor', [0.5 0.7 0.5])
% geoshow('worldrivers.shp','Color', 'blue')
% scatterm([d.lat],[d.lon],200,sum(nb_date,2),'filled'); colorbar;
% 
% figure; imagesc(datenum(datetime('2018-01-01'):datetime('2018-12-31')),1:numel(d),nb_date); 
% yticklabels({d.name}); yticks(1:numel(d)); xlabel('Day of 2018'); datetick('x'); colorbar;
% title('Number of data point per day')


%% Set the regular grid
dt=1/24/4/3;
tnum = datenum(start_date):dt:datenum(end_date);
tdate = datetime(tnum,'convertfrom','datenum');
f_tid = @(x) round((datenum(x)-tnum(1))/dt+1);

for i_d=1:numel(d)
    id = f_tid(d(i_d).time);
    
    for q=1:numel(quantity)
        tmp = nan(numel(tnum),size(d(i_d).DBZH,2));
        tmp(id,:) = d(i_d).(quantity{q});
        d(i_d).(quantity{q})=tmp;
    end
    d(i_d).time = tdate;
end

%% 3. Load Sunrise
% % Downalod sunrise/sunset time at radar location
% time = start_date-1:end_date+1;
% data={};
% for i_d=1:numel(d)
%     data(i_d).lat = d(i_d).lat;
%     data(i_d).lon = d(i_d).lon;
%     data(i_d).name = d(i_d).name;
%     for i_t=1:numel(time)
%         dd = webread(['https://api.sunrise-sunset.org/json?lat=' num2str(d(i_d).lat) '&lng=' num2str(d(i_d).lon) '&date=' datestr(time(i_t),'yyyy-mm-dd')]);
%         assert(strcmp(dd.status,'OK'))
%         data(i_d).sunrs{i_t} = dd.results;
%         % pause(5)
%     end
% end
% save('data/sunrisesunset.mat','time','data')

% Match sunset/sunrise time of radars with d struct
load('data/sunrisesunset.mat')
time2=datestr(time,'yyyy-mmm-dd');
sunr = strings(numel(d), numel(time));
suns = strings(numel(d), numel(time));
for i_d=1:numel(d)
    id = find(strcmp({data.name},d(i_d).name));
    assert(id~=0,'no data for this radar')

    for i_t=1:numel(time)
        sunr(i_d,i_t) = [time2(i_t,:) ' ' data(i_d).sunrs{i_t}.sunrise];%civil_twilight_begin];
        suns(i_d,i_t) = [time2(i_t,:) ' ' data(i_d).sunrs{i_t}.sunset];%civil_twilight_end];
    end
end
sunr = datetime(sunr);
suns = datetime(suns);
for i_d=1:numel(d)
    i=find(strcmp(d(i_d).name,{data.name}));
    [~,Locb] = ismember(dateshift(d(i_d).time,'start', 'day','nearest'),time-0.5);
    d(i_d).sunset = suns(i_d,Locb-1)';
    d(i_d).sunrise = sunr(i_d,Locb)';
end

% Check that it works fine
% figure; hold on;
% plot(d(i_d).time,d(i_d).time,'-r')
% plot(d(i_d).time,d(i_d).sunset,'--k')
% plot(d(i_d).time,d(i_d).sunrise,'--k')


% Apply day mask // can also remove data if needed (later?)
fn=fieldnames(d);
for i_d=1:numel(d)
    d(i_d).day = d(i_d).time'<d(i_d).sunset | d(i_d).time'>d(i_d).sunrise;
%     for i_fn=1:numel(fn)
%         if length(d(i_d).(fn{i_fn}))==length(t_nan)
%             d(i_d).(fn{i_fn})(t_nan,:)=[];
%         end
%     end
end

% Observed data
% CleaningV(d)

%% Automatic cleaning of rain: NOT USED!!!
% Necessery for cleaning raing
for i_d=1:numel(d)
    lim = C.FirstBinOkInDBHZ(strcmp(C.Name,d(i_d).name));
    d(i_d).DBZH(:,1:lim-1)=NaN;
    d(i_d).rain = false(size(d(i_d).DBZH,1),1);
end

% TEST1: Clean Rain with exteranl Rain Data
% 
% Load rain from CDS
rain.file='C:\Users\rnussba1\switchdrive\WeatherRadar\ECMWF\cds.services.cdm_translate-1558601430.4342973-7636-5-be9a0ae2-847b-4096-a3da-5f83f6de5bad_ERA5_REANALYSIS_tcrw_NON_CDM.nc';
% rain.info=ncinfo(rain.file);
rain.tcrw = ncread(rain.file,'tcrw_NON_CDM');
rain.lat=ncread(rain.file,'lat');
rain.lon=ncread(rain.file,'lon');
[rain.LAT,rain.LON] = meshgrid(rain.lat, rain.lon);
rain.time = datetime(ncread(rain.file,'time'), 'ConvertFrom', 'posixtime');

% figure; hold on;
% worldmap([42.5 55.5], [-6 15.5]);  
% geoshow( 'landareas.shp', 'FaceColor', [0.5 0.7 0.5])
% geoshow('worldrivers.shp','Color', 'blue')
% circlem([d.lat],[d.lon],[d.maxrange],'facecolor','r')
% surfm(LAT+.25/2,LON+.25/2,zeros(size(LAT)),'FaceColor','none','EdgeColor','k','LineStyle','-');

% use a subgrid to figure out how much much rain there is over the scanning
% range of each radar. (weighted avg per area)
lat2=linspace(rain.lat(1),rain.lat(end),1000);
lon2=linspace(rain.lon(1),rain.lon(end),1000);
[LAT2,LON2] = meshgrid(lat2, lon2);
ID = reshape(1:numel(rain.LAT),size(rain.LAT));
ID2 = interp2(rain.LAT,rain.LON,ID,LAT2,LON2,'nearest'); 

% Comupte avg rain for each radar
for i_d=1:numel(d)
    dist = lldistkm([d(i_d).lat d(i_d).lon] ,[LAT2(:) LON2(:)]);
    F = reshape(histcounts(ID2(dist(:)<d(i_d).maxrange),.5:(ID(end)+.5)),numel(rain.lon),numel(rain.lat));
    F = F ./ sum(sum(F,1),2);

    tmp = reshape(sum(sum( repmat(F,1,1,size(rain.tcrw,3)) .* rain.tcrw ,1),2),1,[]);
    d(i_d).tcrw = interp1(rain.time,tmp',d(i_d).time,'nearest','extrap');
    
    % rainCS = tmp > 0.0029; % this value is calibrated for bejab
    % rainCS = movmax(rainCS,5);
    % d(i_d).rain = logical(interp1(rain.time,double(rainCS'),d(i_d).time,'nearest','extrap'));
end
clear ID ID2 lon2 lat2 LAT2 LON2 tmp F

% 
% % load('data/bejab_rain.mat');
% % j=0:.0001:.1;
% % for i=1:numel(j)
% %     rainCS = tmp > j(i);
% %     rainCS = movmax(rainCS,3);
% %     rain = logical(interp1(time,double(rainCS'),d(i_d).time,'nearest','extrap'));
% %     r(i)=mean(bejab_rain==rain);
% % end
% % figure; plot(r)
% % [~,i]=max(r);
% % j(i)


% TEST 2: Clearning rain with threashold

% Manual editing for a single radar (bejab)
% Test different quantile value threshold. and look at the % of rain
% classify with this threashold only. median (quantile 50) is a good one
% (40-50% classify with a threashold of ~0. 
% load('data/bejab_rain.mat');
% i_d=1;
% figure;
% a=sort(d(i_d).DBZH,2);
% for q=1:16
%     b=a(:,q);
%     %th=quantile(a(~bejab_rain),.95);
%     subplot(4,4,q); hold on; 
%     histogram(b(bejab_rain & ~d(i_d).rain)); histogram(b(~bejab_rain & ~d(i_d).rain))
%     %plot(th,0,'xr')
%     %legend(num2str(mean(a(bejab_rain)>th)*100))
% end

for i_d=1:numel(d)
    a=sort(d(i_d).DBZH,2,'descend','MissingPlacement','last');
    d(i_d).rain = false(size(d(i_d).DBZH,1),1);
    % DBZH>0 db for the 5th highest value 

    d(i_d).rain(a(:,5)>10)=true;

    % add padding of 30min
    d(i_d).rain = movmax(d(i_d).rain,13);
    
%     tmp = d(i_d).DBZH;
%     tmp(d(i_d).rain,:)=nan;
%     tmp2=tmp-repmat(nanmean(tmp),numel(d(i_d).time),1);
%     lim = C.FirstBinOkInDBHZ(strcmp(C.Name,d(i_d).name));
%     d(i_d).day = -10>(nanmean(tmp2(:,lim:(end-lim)/2+lim),2) - nanmean(tmp2(:,(end-lim)/2+lim+1:end),2));
end

%% Scattering
% Determine which bin are contaminated by scatter with the mean
% MDBZH=nan(25,numel(d));
% for i_d=1:numel(d)
%     MDBZH(1:size(d(i_d).DBZH,2),i_d)=nanmean(d(i_d).DBZH);
% end
% % Visual inspection
% figure; hold on; imagesc(MDBZH);
% xticklabels({d.name}); xticks(1:numel(d)); xtickangle(90)
% plot([d.scatter_lim],'xk')
% 
% This is now loaded with cleaning.xlsx and dealt with in the next section

%% Cleaning manually everything
% Copy all line below. Intrusive method -> create a new variable?
for i_d=1:numel(d)
    d(i_d).dens2 = d(i_d).eta;
    d(i_d).dens2(:,1:d(i_d).scatter_lim) = repmat(d(i_d).eta(:,d(i_d).scatter_lim),1,d(i_d).scatter_lim);  
    d(i_d).dens2 = d(i_d).dens2/11;
    d(i_d).dens2(d(i_d).day,:) = NaN;
end

% Cleaning manually
dc = CleaningV(d);

% save as ('data/dc_corr');

%% Automatic post-cleaning
load('data/dc_corr');

C=readtable('data/Cleaning.xlsx');

for i_d=1:numel(dc)
    
    dc(i_d).alt=dc(i_d).interval*(1:dc(i_d).levels)-dc(i_d).interval/2;
    
    % Expend to 5000m
    dc(i_d).dens3 = [ dc(i_d).dens2 nan(size(dc(i_d).dens2,1),25-size(dc(i_d).dens2,2))];

    % Remove noise level element
    i_dd = strcmp(C.Name,dc(i_d).name);
    tmp = datetime(C.NoiseDate{i_dd},'format','dd/MM/yyyy');
    noise_range = tmp-.5<dc(i_d).time & dc(i_d).time < tmp+.5;
    noise = dc(i_d).dens3(noise_range,:);
    nanmean(noise)
    
    dc(i_d).dens3(any(isnan(dc(i_d).dens2(:, dc(i_d).scatter_lim:dc(i_d).scatter_lim+5)),2),:)=NaN;

    id=~all(isnan(dc(i_d).dens2),2);
    dc(i_d).dens3(id,25)=0;
    dc(i_d).dens3=fillmissing(dc(i_d).dens3,'pchip',2);
end




for i_d=numel(d):-1:1
    id = strcmp(C.Name,d(i_d).name);
    if ~C.Keep(id)
        d(i_d)=[];
    else
        d(i_d).scatter_lim=C.FirstBinOkInEta(id);
        lim = C.FirstBinOkInDBHZ(strcmp(C.Name,d(i_d).name));
        id = strcmp(C.Name,d(i_d).name);
        i = C.FirstBinOk_copyItBelowForScattering(id);
    end
end

figure; hold on; legend
for i_d=1:numel(dc)
    plot(nanmean(dc(i_d).dens3),100:200:4900,'DisplayName',dc(i_d).name)
end

%dc=rmfield('levels','interval');













%% Get elevatation around

files={'gt30w020n90'};%,'gt30e020n90'}; % ,'gt30e020n40','gt30w020n40'
[DEM,GeoCellRef] = geotiffread(['C:\Users\rnussba1\Documents\MATLAB\bin\' files{1}]);
GeoCellRef_lon = linspace(GeoCellRef.LongitudeLimits(1),GeoCellRef.LongitudeLimits(2),GeoCellRef.RasterSize(2));
GeoCellRef_lat = linspace(GeoCellRef.LatitudeLimits(2),GeoCellRef.LatitudeLimits(1),GeoCellRef.RasterSize(1));
[GeoCellRef_LON,GeoCellRef_LAT] = meshgrid(GeoCellRef_lon,GeoCellRef_lat);

figure('position',[0 0 800 600]);
worldmap([floor(min([d.lat])) ceil(max([d.lat]))], [floor(min([d.lon])) ceil(max([d.lon]))]);
geoshow(DEM,GeoCellRef,'DisplayType','texturemap'); demcmap(DEM)

% Set water level to 0
DEM(DEM<100)=0;

% projection angle
alpha=-90:10:90;

% DLAT, DLON
dlat = lldistkm([d.lat; d.lon]', [d.lat; d.lon]' + repmat([GeoCellRef.CellExtentInLatitude 0],numel(d),1));
dlon = lldistkm([d.lat; d.lon]', [d.lat; d.lon]' + repmat([0 GeoCellRef.CellExtentInLongitude],numel(d),1));

% volume in km^3
V = nan(numel(d),numel(alpha));

for i_d=1:numel(d)
    % pre-select coordinate within 1° lat and lon to reduce computational
    % time of lldistkm.
    id=false(GeoCellRef.RasterSize);
    id(abs(GeoCellRef_lat-d(i_d).lat)<1 , abs(GeoCellRef_lon-d(i_d).lon)<1) = true;
    distkm = lldistkm([d(i_d).lat d(i_d).lon],[GeoCellRef_LAT(id) GeoCellRef_LON(id)]);
    distkm_id = distkm<d(i_d).maxrange;
    id_1 = find(id);
    id(id_1(~distkm_id))=0;
    
    % Find DEM of radar -> ELEVATION INCORRECT COMPARED TO GOOGLE MAP
    % ~0-20m (more for higher altitude)
%     [~,i_min] = min(distkm); 
%     d(i_d).heightDEM = DEM(id_1(i_min));

    % We imrotate only the matrix which contain value within distance.
    % we subsample to the rectangle encompasing the radar scan
    id_any_1 = any(id,2); id_any_2 = any(id);
    id_s = id(id_any_1, id_any_2);
    DEM_id = DEM(id_any_1, id_any_2);
    DEM_id(~id_s)=0;
    
    for i_alpha=1:numel(alpha)
        J = imrotate(id_s,alpha(i_alpha));
        DEM_id_r = imrotate(DEM_id,alpha(i_alpha));
        V(i_d,i_alpha) = sum( (d(i_d).height-max(double(DEM_id_r)))/1000 .* sum(J)) *dlat(i_d)*dlon(i_d);
%         subplot(5,4,i_alpha); 
%         % imagesc(DEM_id_r); 
%         plot((d(i_d).height-max(double(DEM_id_r))))
%         title(['a=' num2str(alpha(i_alpha)) '   V=' num2str(V(i_d,i_alpha))])
    end
    d(i_d).VolBelow=V(i_d,:);
%     
    figure; hold on; title(d(i_d).name)
    %scatter3(GeoCellRef_LAT(:), GeoCellRef_LON(:), DEM(:),'.k')
    surf(GeoCellRef_LON(id_any_1, id_any_2), GeoCellRef_LAT(id_any_1, id_any_2), DEM_id,'EdgeColor','none')
    
    plot3(d(i_d).lon, d(i_d).lat, d(i_d).height,'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10);
    xlabel('lat'); ylabel('lon'); zlabel('elev');
    view(3)
end

V2=V./repmat(pi*[d.maxrange]'.^2*0.2,1,numel(alpha));

% figure;
figure; 
for i_d=1:numel(d)
    subplot(6,7,i_d); 
    x1=deg2rad(alpha); y1=V2(i_d,:); x1(V2(i_d,:)<=0)=NaN; y1(V2(i_d,:)<=0)=NaN;
    polarplot(x1,abs(y1))
    hold on;
    x1=deg2rad(alpha); y1=V2(i_d,:); x1(V2(i_d,:)>=0)=NaN; y1(V2(i_d,:)>=0)=NaN;
    polarplot(x1,abs(y1));
    ax=gca;
    ax.ThetaDir = 'clockwise';
    ax.ThetaLim=[-90 90];
    ax.ThetaZeroLocation='top';
    rlim([0 1.2])
    %thetaticks('');
    title(d(i_d).name)
end

% Find elevation with API OPEN_ELEVATION. NOT ALWAYS CORRECT. EXECEL
% spreadsheet has corrected value with https://www.freemaptools.com/elevation-finder.htm
strloc = join(strcat( strcat(cellstr(num2str([d.lat]')),',') , cellstr(num2str([d.lon]')) )','|');
strloc = replace(strloc{1},' ','');
tmp = webread(['https://api.open-elevation.com/api/v1/lookup?locations=' strloc ]); % num2str(d(i_d).lat) ',' num2str(d(i_d).lon)
tmp = jsondecode('{"results": [{"latitude": 51.1917, "elevation": -1, "longitude": 3.0642}, {"latitude": 49.9143, "elevation": 546, "longitude": 5.5056}, {"latitude": 49.6583, "elevation": 861, "longitude": 13.8178}, {"latitude": 49.5011, "elevation": 744, "longitude": 16.7885}, {"latitude": 53.564, "elevation": 1, "longitude": 6.7483}, {"latitude": 54.0044, "elevation": 93, "longitude": 10.0468}, {"latitude": 51.1246, "elevation": 225, "longitude": 13.7686}, {"latitude": 49.5407, "elevation": 774, "longitude": 12.4028}, {"latitude": 53.3394, "elevation": 4, "longitude": 7.025}, {"latitude": 51.4055, "elevation": 153, "longitude": 6.9669}, {"latitude": 51.3112, "elevation": 560, "longitude": 8.802}, {"latitude": 47.8736, "elevation": 1487, "longitude": 8.0036}, {"latitude": 52.4601, "elevation": 58, "longitude": 9.6945}, {"latitude": 48.0431, "elevation": 668, "longitude": 10.2204}, {"latitude": 50.5001, "elevation": 844, "longitude": 11.135}, {"latitude": 50.1097, "elevation": 553, "longitude": 6.5485}, {"latitude": 49.9859, "elevation": 204, "longitude": 8.714}, {"latitude": 52.6486, "elevation": 145, "longitude": 13.8578}, {"latitude": 54.1757, "elevation": 2, "longitude": 12.0581}, {"latitude": 48.1747, "elevation": 633, "longitude": 12.1018}, {"latitude": 48.5853, "elevation": 733, "longitude": 9.7828}, {"latitude": 52.1601, "elevation": 153, "longitude": 11.1761}, {"latitude": 55.1127, "elevation": 162, "longitude": 14.8875}, {"latitude": 55.1731, "elevation": 2, "longitude": 8.552}, {"latitude": 57.4893, "elevation": 85, "longitude": 10.1365}, {"latitude": 55.3262, "elevation": 31, "longitude": 12.4493}, {"latitude": 56.024, "elevation": 120, "longitude": 10.0246}, {"latitude": 59.3977, "elevation": 33, "longitude": 24.6021}, {"latitude": 58.4823, "elevation": 131, "longitude": 25.5187}, {"latitude": 36.8325, "elevation": 479, "longitude": -2.08222}, {"latitude": 39.4289, "elevation": 605, "longitude": -6.28528}, {"latitude": 41.4081, "elevation": 613, "longitude": 1.88472}, {"latitude": 43.1689, "elevation": 592, "longitude": -8.52694}, {"latitude": 41.9956, "elevation": 873, "longitude": -4.60278}, {"latitude": 28.0186, "elevation": 1760, "longitude": -15.6144}, {"latitude": 40.1758, "elevation": 702, "longitude": -3.71361}, {"latitude": 36.6133, "elevation": 1131, "longitude": -4.65917}, {"latitude": 38.2644, "elevation": 1187, "longitude": -1.18972}, {"latitude": 39.3797, "elevation": 109, "longitude": 2.785}, {"latitude": 43.4625, "elevation": 867, "longitude": -6.30194}, {"latitude": 37.6875, "elevation": 501, "longitude": -6.33444}, {"latitude": 43.4033, "elevation": 554, "longitude": -2.84194}, {"latitude": 39.1761, "elevation": 166, "longitude": -0.25211}, {"latitude": 41.7339, "elevation": 811, "longitude": -0.54583}, {"latitude": 60.9039, "elevation": 181, "longitude": 27.1081}, {"latitude": 61.7673, "elevation": 157, "longitude": 23.0764}, {"latitude": 61.907, "elevation": 186, "longitude": 29.7977}, {"latitude": 60.1285, "elevation": 113, "longitude": 21.6434}, {"latitude": 62.8626, "elevation": 107, "longitude": 27.3815}, {"latitude": 67.1391, "elevation": 0, "longitude": 26.8969}, {"latitude": 64.7749, "elevation": 93, "longitude": 26.3189}, {"latitude": 60.2706, "elevation": 160, "longitude": 24.869}, {"latitude": 63.1048, "elevation": 0, "longitude": 23.8209}, {"latitude": 50.1358, "elevation": 69, "longitude": 1.83472}, {"latitude": 42.1297, "elevation": 40, "longitude": 9.49639}, {"latitude": 50.1283, "elevation": 189, "longitude": 3.81194}, {"latitude": 47.3553, "elevation": 588, "longitude": 4.77583}, {"latitude": 44.3231, "elevation": 293, "longitude": 4.76222}, {"latitude": 44.8314, "elevation": 46, "longitude": -0.69167}, {"latitude": 47.0586, "elevation": 161, "longitude": 2.35944}, {"latitude": 48.9272, "elevation": 152, "longitude": -0.14944}, {"latitude": 46.6986, "elevation": 157, "longitude": 0.06556}, {"latitude": 43.2167, "elevation": 616, "longitude": 6.37278}, {"latitude": 45.1044, "elevation": 328, "longitude": 1.36944}, {"latitude": 45.29, "elevation": 1115, "longitude": 3.70944}, {"latitude": 43.9906, "elevation": 664, "longitude": 2.60972}, {"latitude": 43.6247, "elevation": 127, "longitude": -0.60917}, {"latitude": 47.3686, "elevation": 904, "longitude": 7.01917}, {"latitude": 48.7158, "elevation": 282, "longitude": 6.58167}, {"latitude": 43.8061, "elevation": 64, "longitude": 4.50278}, {"latitude": 46.0678, "elevation": 899, "longitude": 4.44528}, {"latitude": 42.9183, "elevation": 688, "longitude": 2.865}, {"latitude": 48.4608, "elevation": 96, "longitude": -4.43}, {"latitude": 43.5744, "elevation": 163, "longitude": 1.37611}, {"latitude": 48.7739, "elevation": 168, "longitude": 2.0075}, {"latitude": 47.3375, "elevation": 66, "longitude": -1.65639}, {"latitude": 48.4622, "elevation": 151, "longitude": 4.30944}, {"latitude": 45.8835, "elevation": 246, "longitude": 17.2009}, {"latitude": 45.5027, "elevation": 87, "longitude": 18.5613}, {"latitude": 52.9528, "elevation": 6, "longitude": 4.79061}, {"latitude": 51.8369, "elevation": 2, "longitude": 5.1381}, {"latitude": 50.3942, "elevation": 393, "longitude": 20.0797}, {"latitude": 54.3843, "elevation": 135, "longitude": 18.4563}, {"latitude": 52.4052, "elevation": 92, "longitude": 20.9609}, {"latitude": 50.892, "elevation": 662, "longitude": 16.0395}, {"latitude": 52.4133, "elevation": 94, "longitude": 16.7971}, {"latitude": 50.1517, "elevation": 322, "longitude": 18.7267}, {"latitude": 50.1141, "elevation": 207, "longitude": 22.037}, {"latitude": 53.7903, "elevation": 113, "longitude": 15.8311}, {"latitude": 56.3675, "elevation": 190, "longitude": 12.8517}, {"latitude": 59.6544, "elevation": 50, "longitude": 17.9463}, {"latitude": 57.3034, "elevation": 58, "longitude": 18.4003}, {"latitude": 61.5771, "elevation": 0, "longitude": 16.7144}, {"latitude": 67.7088, "elevation": 0, "longitude": 20.6178}, {"latitude": 56.2955, "elevation": 106, "longitude": 15.6103}, {"latitude": 60.723, "elevation": 0, "longitude": 14.8776}, {"latitude": 65.4309, "elevation": 0, "longitude": 21.865}, {"latitude": 63.6395, "elevation": 64, "longitude": 18.4019}, {"latitude": 63.295, "elevation": 251, "longitude": 14.7591}, {"latitude": 58.2556, "elevation": 136, "longitude": 12.826}, {"latitude": 58.1059, "elevation": 206, "longitude": 15.9363}, {"latitude": 46.0678, "elevation": 937, "longitude": 15.2849}, {"latitude": 46.098, "elevation": 1004, "longitude": 14.2282}, {"latitude": 48.2561, "elevation": 577, "longitude": 17.1531}, {"latitude": 48.7829, "elevation": 1236, "longitude": 20.9873}, {"latitude": 49.2717, "elevation": 1387, "longitude": 19.2493}, {"latitude": 48.2404, "elevation": 613, "longitude": 19.2574}]}');
[tmp.results.elevation]';
%% Export for online visualization
for i_d=1:numel(d)
    fileID = fopen(['data/exportdDens_' d(i_d).name '.json'],'w');
    a=struct();
    a.x=datestr(d(i_d).time,'yyyy-mm-dd HH:MM');
    a.y=d(i_d).interval*(1/2:double(d(i_d).levels));
    b = log10(d(i_d).dens);
    b(d(i_d).dens==0)=0;
    b(b<0)=0;
    a.data=round(b,2)';
    fprintf(fileID,jsonencode(a));
    fclose(fileID);
end
fileID = fopen('data/exportdRadars.json','w');
fprintf(fileID,jsonencode({d.name}));
fclose(fileID);


%% Test DBZH -> dens
i_d=14; % seems like a good case study

% eta-> dens
figure; hold on; 
plot([0 max(d(i_d).dens(:))],[0 max(d(i_d).dens(:))],'--r');
plot(d(i_d).eta(:)/11,d(i_d).dens(:),'.k')

sum(~isnan(d(i_d).eta(isnan(d(i_d).dens(:))))) % NO NaN in dens with value in eta
mean(isnan(d(i_d).eta(d(i_d).dens(:)==0))) % very few dens=0 -> eta=nan
mean(d(i_d).eta(d(i_d).dens(:)==0)==0) % very few dens=0 -> eta=0

id = d(i_d).sd_vvp(:)>2;
figure; 
plot(d(i_d).DBZH(id),d(i_d).sd_vvp(id),'.k')
legend(num2str(sum(id)))

scatter(d(i_d).eta(:)/11,d(i_d).dens(:),[],d(i_d).sd_vvp(:)<2,'.')


% DBZH -> eta
figure; hold on;
a=361*10.^(d(i_d).DBZH(:)/10);
b=d(i_d).eta(:);
plot([0 max(d(i_d).eta(:))],[0 max(d(i_d).eta(:))],'--r');
plot(a,b,'.k');
axis([0 max(d(i_d).eta(:)) 0 max(d(i_d).eta(:))])

% Rain between DBZH and eta
figure; plot(d(i_d).eta(:),d(i_d).DBZH(:),'.')

figure; hold on;
histogram(d(i_d).DBZH(isnan(d(i_d).eta(:))))
histogram(d(i_d).DBZH(~isnan(d(i_d).eta(:))))
