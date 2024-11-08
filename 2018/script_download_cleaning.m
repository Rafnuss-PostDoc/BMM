addpath('functions/')

%% 1. Download value
start_date='01-Jan-2020 00:00:00';
end_date='01-Jan-2021 00:00:00';
quantity = {'dens','ff','dd','DBZH','eta','sd_vvp'};
folder_root = 'G:\Raphael_temp\'; % 'C:\Users\rnussba1\Documents\BMM\data\'

d = download_vp(start_date,end_date, quantity, folder_root);

save('data/d2020all.mat','d','start_date','end_date','quantity','-v7.3');

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

%% Check 2019
start_date='01-Jan-2019 00:00:00';
end_date='01-Jan-2020 00:00:00';
quantity = {'dens','ff','dd','DBZH','eta','sd_vvp'};

d = download_vp(start_date,end_date, quantity);

save('data/d2019all.mat','d','start_date','end_date','quantity','-v7.3');



tmp_t = (datenum(start_date)-.5):(datenum(end_date)+.5);
tmp = nan(numel(tmp_t)-1, numel(d));
for i_d=1:numel(d)
    tmp(:,i_d) = histcounts (round(datenum(d(i_d).time(~all(isnan(d(i_d).dens),2)))),tmp_t);
end
figure;
imagesc(tmp')
datetick('x'); yticks(1:numel(d)); yticklabels({d.name});

%% 2. Delete entire radar 
load('data/d2018all.mat'); addpath('functions/')

% Observed data
% CleaningV(d)

% Load information about cleaning
C=readtable('data/Cleaning.xlsx');
% OPERA_RADARS_DB = struct2table(webread('https://www.eumetnet.eu/wp-content/themes/aeron-child/observations-programme/current-activities/opera/database/OPERA_Database/OPERA_RADARS_DB.json'));
% OPERA_RADARS_DB.Name=OPERA_RADARS_DB.odimcode;
% OPERA_RADARS_DB(~ismember(OPERA_RADARS_DB.Name,C.Name),:)=[];
% C = join(C,OPERA_RADARS_DB,'Keys','Name');
% C.heightofstation=str2double(string(C.heightofstation));
% C.heightantenna=str2double(string(C.heightantenna));

% Remove radars with wrong, imcomplete data: 76 radars overs 123
for i_d=numel(d):-1:1
    id = strcmp(C.Name,d(i_d).name);
    if ~C.Keep(id)
        d(i_d)=[];
    else
        d(i_d).scatter_lim=C.FirstBinOkInEta(id);
        d(i_d).scatter_lim_DBZH=C.FirstBinOkInDBHZ(id);
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
% Median resolution is 5 minutes
% cellfun(@(x) median(diff(x)), {d.time},'UniformOutput',false)

dt=1/24/4/3; % All radar have a 5 minutes delta t. 
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
    
    % altitude axis
    d(i_d).alt=d(i_d).interval*(1:d(i_d).levels)-d(i_d).interval/2;
end

%% 3. Load Sunrise
val = jsondecode(fileread('sunriseset_astral\sunriseset_latlon.json'));
time = datetime(start_date):datetime(end_date);
sunrs=strings(numel(val),numel(val{1}),numel(val{1}{1}));
for i=1:numel(val)
    for j=1:numel(val{i})
        try
            sunrs(i,j,:)=val{i}{j};
        end
    end
end
sunrs = datetime(sunrs);

for i_d=1:numel(d)
    i_d2 = find(strcmp(d(i_d).name, {d.name}));
    [~,Locb] = ismember(dateshift(d(i_d).time,'start', 'day','nearest'),time);
    d(i_d).dawn = sunrs(i_d2,Locb,1)';
    d(i_d).sunrise = sunrs(i_d2,Locb,2)';
    d(i_d).sunset = sunrs(i_d2,max(Locb-1,1),3)';
    d(i_d).dusk = sunrs(i_d2,max(Locb-1,1),4)';
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
% Clean Rain with exteranl Rain Data
% % 
% % Load rain from CDS
% rain.file='C:\Users\rnussba1\switchdrive\WeatherRadar\ECMWF\cds.services.cdm_translate-1558601430.4342973-7636-5-be9a0ae2-847b-4096-a3da-5f83f6de5bad_ERA5_REANALYSIS_tcrw_NON_CDM.nc';
% % rain.info=ncinfo(rain.file);
% rain.tcrw = ncread(rain.file,'tcrw_NON_CDM');
% rain.lat=ncread(rain.file,'lat');
% rain.lon=ncread(rain.file,'lon');
% [rain.LAT,rain.LON] = meshgrid(rain.lat, rain.lon);
% rain.time = datetime(ncread(rain.file,'time'), 'ConvertFrom', 'posixtime');
% 
% % figure; hold on;
% % worldmap([42.5 55.5], [-6 15.5]);  
% % geoshow( 'landareas.shp', 'FaceColor', [0.5 0.7 0.5])
% % geoshow('worldrivers.shp','Color', 'blue')
% % circlem([d.lat],[d.lon],[d.maxrange],'facecolor','r')
% % surfm(LAT+.25/2,LON+.25/2,zeros(size(LAT)),'FaceColor','none','EdgeColor','k','LineStyle','-');
% 
% % use a subgrid to figure out how much much rain there is over the scanning
% % range of each radar. (weighted avg per area)
% lat2=linspace(rain.lat(1),rain.lat(end),1000);
% lon2=linspace(rain.lon(1),rain.lon(end),1000);
% [LAT2,LON2] = meshgrid(lat2, lon2);
% ID = reshape(1:numel(rain.LAT),size(rain.LAT));
% ID2 = interp2(rain.LAT,rain.LON,ID,LAT2,LON2,'nearest'); 
% 
% % Comupte avg rain for each radar
% for i_d=1:numel(d)
%     dist = lldistkm([d(i_d).lat d(i_d).lon] ,[LAT2(:) LON2(:)]);
%     F = reshape(histcounts(ID2(dist(:)<d(i_d).maxrange),.5:(ID(end)+.5)),numel(rain.lon),numel(rain.lat));
%     F = F ./ sum(sum(F,1),2);
% 
%     tmp = reshape(sum(sum( repmat(F,1,1,size(rain.tcrw,3)) .* rain.tcrw ,1),2),1,[]);
%     d(i_d).tcrw = interp1(rain.time,tmp',d(i_d).time,'nearest','extrap');
%     
%     % rainCS = tmp > 0.0029; % this value is calibrated for bejab
%     % rainCS = movmax(rainCS,5);
%     % d(i_d).rain = logical(interp1(rain.time,double(rainCS'),d(i_d).time,'nearest','extrap'));
% end
% clear ID ID2 lon2 lat2 LAT2 LON2 tmp F

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

% deoft has a lot of noise in eta. Take dens and clean it. 
i_d=strcmp({dc.name},'depro'); %deoft, dehnr, depro
dc(i_d).dens2=dc(i_d).dens;
dc(i_d).dens2(all(isnan(dc(i_d).dens3),2),:) = NaN;
for i_t=1:numel(dc(i_d).time)
    firt=find(dc(i_d).dens(i_t,:)>0,1,'last');
    id=dc(i_d).dens2(i_t,1:firt)==0 | isnan(dc(i_d).dens2(i_t,1:firt));
    dc(i_d).dens2(i_t,id)= dc(i_d).dens3(i_t,id);
end

%Save 
for i_d=1:numel(dc)
    dens2{i_d} = dc(i_d).dens2;
end
% First version: limited edit
% save('data/dc_dens2_v1','dens2');
% Second version: ground clutter elimination, more extensif, manual
% clearning
% save('data/dc_dens2_v2','dens2');


%% Automatic post-cleaning
load('data/dc_dens2_v2');
for i_d=1:numel(dc)
    dc(i_d).dens2 = dens2{i_d};
end
%dc=d;

C=readtable('data/Cleaning.xlsx');

% kernel for data alone, if less than 30min of data
B1 = ones(60/5+1,1);
B1(30/5+1)=0;

% kernel for noise removal
B2 = 1/8*ones(3,3);
B2(2,2)=0;

% interpolation grid
[X,Y]=meshgrid(1:25,1:numel(dc(1).time));

for i_d=1:length(dc)
    % altitude axis
    dc(i_d).alt=dc(i_d).interval*(1:dc(i_d).levels)-dc(i_d).interval/2;
    
    % reshape to same size for all radars
    dc(i_d).dens3 = [ dc(i_d).dens2 nan(size(dc(i_d).dens2,1),25-size(dc(i_d).dens2,2))];
    dc(i_d).ff = [ dc(i_d).ff nan(size(dc(i_d).ff,1),25-size(dc(i_d).ff,2))];
    dc(i_d).dd = [ dc(i_d).dd nan(size(dc(i_d).dd,1),25-size(dc(i_d).dd,2))];
    dc(i_d).sd_vvp = [ dc(i_d).sd_vvp nan(size(dc(i_d).sd_vvp,1),25-size(dc(i_d).sd_vvp,2))];
    
    % Remove data if less than 5 lowest layers present
    dc(i_d).dens3(any(isnan(dc(i_d).dens2(:, dc(i_d).scatter_lim:dc(i_d).scatter_lim+5)),2),:)=NaN;
  
    % Remove for alone data 
    id_nan=any(~isnan(dc(i_d).dens3),2);
    tmp = conv2(id_nan,B1,'same');
    dc(i_d).dens3(id_nan & tmp<(30/5-1),:)=nan;

    % Remove noise level element
    ldc = log10(dc(i_d).dens3+1);
    id_nan = isnan(ldc);
    tmp2=conv2(~id_nan,B2,'same');
    ldc(id_nan)=0;
    idx = ldc > 2*conv2(ldc,B2,'same')./tmp2;
    
    % Interpolate up to 5000m. 
    id=~all(isnan(dc(i_d).dens3),2);
    dc(i_d).dens3(id,25)=0;
    idx2 = repmat(id,1,25) & id_nan;

    idp = ~isnan(dc(i_d).dens3) & ~idx;
    F = scatteredInterpolant(X(idp), Y(idp), log10(dc(i_d).dens3(idp)+1),'natural');
    
    dc(i_d).dens3(idx) = 10.^F(X(idx), Y(idx))-1;
    dc(i_d).dens3(idx2) = 10.^F(X(idx2), Y(idx2))-1;
    dc(i_d).dens3(dc(i_d).dens3<0) = 0;
end

%Save 
for i_d=1:numel(dc)
    dens3{i_d} = dc(i_d).dens3;
end
save('data/dc_dens3','dens3');


save('data/dc_corr','dc','start_date','end_date','quantity','-v7.3')















%% Automatic clearning: Calibration on manually cleaned one
load('data/dc_dens2_v2')
for i_d=1:numel(d)
    d(i_d).dens2 = dens2{i_d};
end


% Calibrate rain detection paramter
parm_rain=nan(5,numel(d));
for i_d=1:numel(d)
    fun_cal = @(x) autoCleanRain_calibration(d(i_d),x);
    
%     lb = [1 -10 1 -10 1];
%     ub = [8 20 8 20 30];
%     x0 = [6 10 4 0 12];
%     IntCon = 1:5;
%     options = optimoptions('ga','Display','iter','PlotFcn',@gaplotbestf,'FunctionTolerance',1e-1);
%     parm_rain(:,i_d)  = ga(fun_cal,5,[],[],[],[],lb,ub,[],IntCon,options);

    d(i_d).rain = autoCleanRain(d(i_d),parm_rain(1,i_d),parm_rain(2,i_d),parm_rain(3,i_d),parm_rain(4,i_d),parm_rain(5,i_d));
    
    % Convert eta to density and apply mask
    d(i_d).dens3 = d(i_d).eta/11;
    d(i_d).dens3(d(i_d).day | d(i_d).rain,:) = nan;

    % reshape to same size for all radars
    d(i_d).dens3 = [ d(i_d).dens3 nan(size(d(i_d).dens3,1),25-size(d(i_d).dens3,2))];

end



% Calibrate interpolat detection paramter
parm_interp=nan(2,numel(d));
for i_d=1:numel(d)

    fun_cal = @(x) autoCleanInterp_calibration(d(i_d),x);
    lb = [1 -10 1 -10 1];
    ub = [8 20 8 20 30];
    x0 = [6 10 4 0 12];
    IntCon = 1:5;
    options = optimoptions('ga','Display','iter','PlotFcn',@gaplotbestf,'FunctionTolerance',1e-1);
    parm_interp(:,i_d)  = ga(fun_cal,5,[],[],[],[],lb,ub,[],IntCon,options);
    d(i_d).interp = autoCleanInterp(d(i_d),parm_interp(1,i_d),parm_interp(2,i_d));
end


save('data/autoCleanX','parm_interp','parm_rain')

% Perform the autoclearning + interpolation

for i_d=1:numel(d)
    % Convert eta to density and apply mask
    d(i_d).dens3 = d(i_d).eta/11;
    d(i_d).dens3(d(i_d).day | d(i_d).rain,:) = nan;

    % reshape to same size for all radars
    d(i_d).dens3 = [ d(i_d).dens3 nan(size(d(i_d).dens3,1),25-size(d(i_d).dens3,2))];
    % di.ff = [ di.ff nan(size(di.ff,1),25-size(di.ff,2))];
    % di.dd = [ di.dd nan(size(di.dd,1),25-size(di.dd,2))];

    % Build the scatter interpolant
    r = 200; % ratio of correlation between vertical and time: 1 => 200m = 5minutes
    [X,Y]=meshgrid(linspace(1,r*25,25),1:numel(di.time));
    idp = ~isnan(d(i_d).dens3) & ~id_noise & ~id_raintop;
    F = scatteredInterpolant(X(idp), Y(idp), log10(d(i_d).dens3(idp)+1),'linear');

    d(i_d).dens3(d(i_d).interp) = 10.^F(X(d(i_d).interp), Y(d(i_d).interp))-1;
    d(i_d).dens3(d(i_d).dens3<0) = 0;
end

CleaningV(d)





%% Correction for flight speed and sd_vvp

r = 200; % ratio of correlation between vertical and time: 1 => 200m = 5minutes
[X,Y]=meshgrid(linspace(1,r*25,25),1:numel(dc(i_d).time));
  
for i_d=1:numel(dc)
    % Threashold speed : NOT APPLY HERE
    % dc(i_d).ff(dc(i_d).ff>25)=nan;
    
    % Compute the n/s and e/w componenent of the flight speed
    dc(i_d).u = dc(i_d).ff .* sind(dc(i_d).dd); % m/s | 0� is north and 90� is west. -> u is east (+) - west (-)
    dc(i_d).v = dc(i_d).ff .* cosd(dc(i_d).dd); % m/s  -> v is north (+) - south (-)

    % Remove data not cleaned for density (e.g. rain). 
    id = isnan(dc(i_d).dens3) | dc(i_d).dens3==0;
    dc(i_d).u(id)=nan;
    dc(i_d).v(id)=nan;
    dc(i_d).sd_vvp2 = dc(i_d).sd_vvp;
    dc(i_d).sd_vvp2(id)=nan;
    
    Fu = scatteredInterpolant(X(~isnan(dc(i_d).u)), Y(~isnan(dc(i_d).u)), dc(i_d).u(~isnan(dc(i_d).u)),'natural','none');
    Fv = scatteredInterpolant(X(~isnan(dc(i_d).v)), Y(~isnan(dc(i_d).v)), dc(i_d).v(~isnan(dc(i_d).v)),'natural','none');
    Fsd_vvp = scatteredInterpolant(X(~isnan(dc(i_d).sd_vvp)), Y(~isnan(dc(i_d).sd_vvp)), dc(i_d).sd_vvp(~isnan(dc(i_d).sd_vvp)),'natural','none');
    % Keep data only if nearby 2 data point are whithin +/-2hr or +/-400m
    tmp=conv2(~isnan(dc(i_d).u),ones(4*4*2+1,2*2+1),'same');
    idp = tmp>2 & isnan(dc(i_d).u) & ~isnan(dc(i_d).dens3);
    dc(i_d).u(idp) = Fu(X(idp), Y(idp));
    tmp=conv2(~isnan(dc(i_d).v),ones(4*4*2+1,2*2+1),'same');
    idp = tmp>2 & isnan(dc(i_d).v) & ~isnan(dc(i_d).dens3);
    dc(i_d).v(idp) = Fv(X(idp), Y(idp));
    tmp=conv2(~isnan(dc(i_d).sd_vvp2),ones(4*4*2+1,2*2+1),'same');
    idp = tmp>2 & isnan(dc(i_d).sd_vvp2) & ~isnan(dc(i_d).dens3);
    dc(i_d).sd_vvp2(idp) = Fsd_vvp(X(idp), Y(idp));

    % OLD INTERPOLATION: CAN DELETE
%     dc(i_d).u = fillmissing(dc(i_d).u,'linear');
%     dc(i_d).v = fillmissing(dc(i_d).v,'linear');
%     dc(i_d).u(tmp<4)=nan;
%     dc(i_d).v(tmp<4)=nan;
    
% MABYE NEEDED FOR THE INTERPOLATION
    % check that at direction is representatife of at least 50% of the bird
    nb_bird = dc(i_d).dens3.* repmat(max(min(dc(1).alt-dc(i_d).heightDEM,200),0)/1000,numel(dc(i_d).time),1);
    nb_bird_2=nb_bird;
    nb_bird_2(isnan(dc(i_d).u))=nan;
    id = nansum(nb_bird_2,2)<.5*nansum(nb_bird,2);
    dc(i_d).u(id,:)=nan;
    dc(i_d).v(id,:)=nan;
    dc(i_d).sd_vvp2(id,:)=nan;
end












%% View result and save
CleaningV(dc)

save('data/dc_corr','dc','start_date','end_date','quantity','-v7.3')

