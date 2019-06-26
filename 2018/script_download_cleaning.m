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
        d(i_d).heightDEM=C.HeightDEM(id);
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

%Save 
for i_d=1:numel(dc)
    dens2{i_d} = dc(i_d).dens2;
end
% First version: limited edit
% save('data/dc_dens2_v2','dens2');
% Second version: ground clutter elimination, more extensif
% save('data/dc_dens2_v2','dens2');


%% Automatic post-cleaning
load('data/dc_dens2_v2');
for i_d=1:numel(d)
    d(i_d).dens2 = dens2{i_d};
end
dc=d;

C=readtable('data/Cleaning.xlsx');
B = 1/8*ones(3,3);
B(2,2)=0;
[X,Y]=meshgrid(1:25,1:numel(dc(1).time));

for i_d=1:length(dc)
    
    dc(i_d).alt=dc(i_d).interval*(1:dc(i_d).levels)-dc(i_d).interval/2;
    
    % Expend to 5000m
    dc(i_d).dens3 = [ dc(i_d).dens2 nan(size(dc(i_d).dens2,1),25-size(dc(i_d).dens2,2))];
    
    % Remove data if less than 5 lowest layers present
    dc(i_d).dens3(any(isnan(dc(i_d).dens2(:, dc(i_d).scatter_lim:dc(i_d).scatter_lim+5)),2),:)=NaN;
    
    % Remove noise level element
%     i_dd = strcmp(C.Name,dc(i_d).name);
%     if ~strcmp(C.NoiseDate{i_dd},'NaN')
%         tmp = datetime(C.NoiseDate{i_dd},'format','dd/MM/yyyy');
%         noise_range = tmp-.5<dc(i_d).time & dc(i_d).time < tmp+.5;
%         noise = dc(i_d).dens3(noise_range,:);
%         dc(i_d).dens3 = dc(i_d).dens3-repmat(nanmean(noise),size(dc(i_d).dens3,1),1);
%         dc(i_d).dens3(dc(i_d).dens3<0)=0;
%     end

    
    ldc = log10(dc(i_d).dens3+1);
    id_nan = isnan(ldc);
    tmp2=conv2(~id_nan,B,'same');
    ldc(id_nan)=0;
    idx = ldc > 2*conv2(ldc,B,'same')./tmp2;
    

    % Interpolate up to 5000m. 
    id=~all(isnan(dc(i_d).dens2),2);
    dc(i_d).dens3(id,25)=0;
    idx2 = repmat(id,1,25) & id_nan;

    idp = ~isnan(dc(i_d).dens3) & ~idx;
    F = scatteredInterpolant(X(idp), Y(idp), log10(dc(i_d).dens3(idp)+1),'natural');
    
    dc(i_d).dens3(idx) = 10.^F(X(idx), Y(idx))-1;
    dc(i_d).dens3(idx2) = 10.^F(X(idx2), Y(idx2))-1;
end

%Save 
for i_d=1:numel(dc)
    dens3{i_d} = dc(i_d).dens3;
end
% save('data/dc_dens3','dens3');


% save('data/dc_corr','dc','start_date','end_date','quantity','-v7.3')




