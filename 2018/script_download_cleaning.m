
%% 1. Download value
start_date='01-Jan-2018 00:00:00';
end_date='01-Jan-2019 00:00:00';
quantity = {'dens','ff','dd','n_all','DBZH','w'};

d = download_vp(start_date,end_date, quantity);

save('data/d2018all.mat','d','start_date','end_date','quantity');

%% 2. Delete entire radar 
load('data/d2018_all.mat'); addpath('functions/')

% Observed data
CleaningV(d)

% Remove radars with wrong, imcomplete data: 76 radars overs 123
radar_rm={'bezav','bg*','ct*','cz*','deemd','defbg','deflg','deisn','dk*','ee*','es*','fi*','frale','frbol','fropo','hr*','nldbl','pl*','pt*','searl','sease','sehud','sekir','sekkr','selul','seosd','seosu','seovi','sevar','sevil','si*','sk*'};

d_rm=[];
for i_d=1:numel(d)
    res = regexp(d(i_d).name,regexptranslate('wildcard',radar_rm));
    if any(cell2mat(res')==1)
        d_rm=[d_rm;i_d];
    end
end
d(d_rm)=[];

nb_date=false(numel(d),365);
for i_d=1:numel(d)
    nb_date(i_d,unique(round(datenum(d(i_d).time)))-datenum('1-1-2018'))= 1;
end

figure; hold on;
worldmap([min([d.lat]) max([d.lat])], [min([d.lon]) max([d.lon])]);  
geoshow( 'landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow('worldrivers.shp','Color', 'blue')
scatterm([d.lat],[d.lon],200,sum(nb_date,2),'filled'); caxis([0 365]); colorbar;

figure; imagesc(datenum(datetime('2018-01-01'):datetime('2018-12-31')),1:numel(d),nb_date); yticklabels({d.name}); yticks(1:numel(d)); xlabel('Day of 2018'); datetick('x')


%% 3. Load Sunrise
time = start_date-1:end_date+1;
data={};
for i_d=1:numel(d)
    data(i_d).lat = d(i_d).lat;
    data(i_d).lon = d(i_d).lon;
    data(i_d).name = d(i_d).name;
    for i_t=1:numel(time)
        dd = webread(['https://api.sunrise-sunset.org/json?lat=' num2str(d(i_d).lat) '&lng=' num2str(d(i_d).lon) '&date=' datestr(time(i_t),'yyyy-mm-dd')]);
        assert(strcmp(dd.status,'OK'))
        data(i_d).sunrs{i_t} = dd.results;
        % pause(5)
    end
end
save('data/sunrisesunset.mat','time','data')

load('data/sunrisesunset.mat')
time2=datestr(time,'yyyy-mmm-dd');
sunr = strings(numel(d), numel(time));
suns = strings(numel(d), numel(time));
for i_d=1:numel(d)
    for i_t=1:numel(time)
        sunr(i_d,i_t) = [time2(i_t,:) ' ' data(i_d).sunrs{i_t}.sunrise];%civil_twilight_begin];
        suns(i_d,i_t) = [time2(i_t,:) ' ' data(i_d).sunrs{i_t}.sunset];%civil_twilight_end];
    end
end
sunr = datetime(sunr);
suns = datetime(suns);
for i_d=1:numel(d)
    i=strcmp(d(i_d).name,{data.name});
    [~,Locb] = ismember(dateshift(d(i_d).time,'start', 'day','nearest'),time-0.5);
    d(i_d).sunset = suns(i,Locb-1)';
    d(i_d).sunrise = sunr(i,Locb)';
end

figure; hold on;
plot(d(i_d).time,d(i_d).time,'-r')
plot(d(i_d).time,d(i_d).sunset,'--k')
plot(d(i_d).time,d(i_d).sunrise,'--k')


%% 3. Remove data during day
for i_d=1:numel(d)
    t_avg_i_d = mean([d(i_d).etime d(i_d).stime],2);
    t_nan = t_avg_i_d<d(i_d).sunset | t_avg_i_d>d(i_d).sunrise;
    d(i_d).time(t_nan)=[];
    d(i_d).stime(t_nan)=[];
    d(i_d).etime(t_nan)=[];
    d(i_d).dens(t_nan,:)=[];
    d(i_d).ff(t_nan,:)=[];
    d(i_d).dd(t_nan,:)=[];
    d(i_d).n_all(t_nan,:)=[];
    d(i_d).DBZH(t_nan,:)=[];
    d(i_d).sunrise(t_nan)=[];
    d(i_d).sunset(t_nan)=[];
    d(i_d).w(t_nan,:)=[];
end

% Observed data
CleaningV(d)