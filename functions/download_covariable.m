function [d, co] = download_covariable(d)

% Meteo
% Data are donwloaded form ECMWF with a python script
% (../ECMWF/script.py). 

co.qa = {'u','cc','v','t','crwc'};
co.qs = {'u100','v100','hcc','mcc','lcc','lsm','t2m','v10','u10','tcc','stl1','slor','anor','isor','sdor','tcw','tcrw','tclw','slt'};

% 1 for each radar location
% The lat-lon are generated with the code below
% dpos=90./(1:900);
% dlat = [d.maxrange]./110.574;
% dlon = [d.maxrange]./(111.320*cosd([d.lat]));
% for i=1:numel(dlat)
%     [~,id] = min(abs(dlat(i)-dpos));
%     dlat(i)=dpos(id);
% end
% for i=1:numel(dlon)
%     [~,id] = min(abs(dlon(i)-dpos));
%     dlon(i)=dpos(id);
% end
% coord = [1:numel(d) ; [d.lat]+dlat ; [d.lon]-dlon ; [d.lat]-dlat ; [d.lon]+dlon ; dlat*2 ; dlon*2]';
% pb with 3, 22, 25, 31, 34, 35, 43, 47, 48, 62
%
% for i_d=1:numel(d)
%     file = ['../ECMWF/output_09-2016_' num2str(i_d) '.nc'];
% 
%     % ncdisp(file);
%     
%     % Time is the same for all files
%     if i_d==1
%         t = datetime(single(ncread(file,'time'))/24+datenum('1900-01-01'), 'ConvertFrom', 'datenum');
%     end
%     
%     lat = flip(double(ncread(file,'latitude')));
%     lon = double(ncread(file,'longitude'));
%     
%     for i_q = 1:numel(quantity)
%         q = flip(ncread(file,quantity{i_q}),2);
% 
%         % create the interpolation
%         F = griddedInterpolant({lon,lat,datenum(t)},q);
%         
%         % fun = @(x,y,z) F(x,y,z);
%         % integral3(fun,d(i_d).lon-0.1,d(i_d).lon+0.1,d(i_d).lat-0.1,d(i_d).lat+0.1,datenum(d(i_d).time)-0.1,datenum(d(i_d).time)+0.1)
%         
%         tmp = [F({d(i_d).lon,d(i_d).lat,datenum(d(i_d).time)})];
%         d(i_d).(quantity{i_q}) = tmp(:);
%     end
% end




% For the region "grid": "0.225/0.225", "area": "europe",
alt={'950','900', '850', '800', '750', '700'};
for i_file=1:numel(alt)
    file = ['../ECMWF/output_20160919-20161010_' alt{i_file} '_europe.nc'];
    % ncdisp(file);
    
    co.lat = double(ncread(file,'latitude'));
    co.lon = double(ncread(file,'longitude'));
    co.time = datetime(single(ncread(file,'time'))/24+datenum('1900-01-01'), 'ConvertFrom', 'datenum');
    
    for i_q = 1:numel(co.qa)
        co.(['a' alt{i_file}]).(co.qa{i_q}) = ncread(file,co.qa{i_q});
        co.qal{i_q} = ncreadatt(file,co.qa{i_q},'long_name');
    end
end

file = ['../ECMWF/output_20160919-20161010_srf_europe.nc'];
for i_q = 1:numel(co.qs)
    co.srf.(co.qs{i_q}) = ncread(file,co.qs{i_q});
    co.qsl{i_q} = ncreadatt(file,co.qs{i_q},'long_name');
end

% Get data at the radar locations
for i_file=1:numel(alt)
    for i_q = 1:numel(co.qa)
        F = griddedInterpolant({co.lon,-co.lat,datenum(co.time)},co.(['a' alt{i_file}]).(co.qa{i_q}));
        for i_d=1:numel(d)
            tmp = F({d(i_d).lon,-d(i_d).lat,datenum(d(i_d).time)});
            if any(isnan(tmp(:)))
                keyboard
            end
            d(i_d).(co.qa{i_q})(:,i_file) = tmp(:);
        end
    end
end

for i_q = 1:numel(co.qs)
    F = griddedInterpolant({co.lon,-co.lat,datenum(co.time)},co.srf.(co.qs{i_q}));
    for i_d=1:numel(d)
        tmp = F({d(i_d).lon,-d(i_d).lat,datenum(d(i_d).time)});
        if any(isnan(tmp(:)))
            keyboard
        end
        d(i_d).(co.qs{i_q}) = tmp(:);
    end
end


% Civil Twilight
% SunriseSunset(d, start_date,end_date)
load('./data/sunrisesunset.mat')
time2=datetime(datestr(time,'yyyy-mmm-dd'));
for i_d=1:numel(d)
    d(i_d).sunrise = NaT(1, numel(d(i_d).time));
    d(i_d).sunset = NaT(1, numel(d(i_d).time));
    
    i=strcmp(d(i_d).name,{data.name});
    
    for i_t=2:numel(time)
        id = d(i_d).time>time(i_t)-1 & d(i_d).time<=time(i_t);

        d(i_d).sunset(id) = [datestr(time2(i_t-1),'yyyy-mmm-dd') ' ' data(i).sunrs{i_t-1}.civil_twilight_end ];
        d(i_d).sunrise(id) = [datestr(time2(i_t),'yyyy-mmm-dd') ' ' data(i).sunrs{i_t}.civil_twilight_begin ];
    end
end


co.qs = {co.qs{:} ,'sunrise','sunset'};
co.qsl = {co.qsl{:} ,'sunrise of the nearest night','sunset of the nearest night'};


% figure;
% latlim = [min(co.lat) max(co.lat)];
% lonlim = [min(co.lon) max(co.lon)];
% 
% worldmap('Europe');
% [lat lon] = meshgrat(latlim,lonlim,[numel(co.lat) numel(co.lon) ]);
% pcolorm(lat, lon, co.u(:,:,end)')
% load coastlines; 
% plotm(coastlat, coastlon)


% for i_d=1:numel(d)
% figure(2);clf; hold on;
% subplot(2,1,1);
% plot(d(i_d).time,d(i_d).u)
% subplot(2,1,2);
% plot(d(i_d).time,d(i_d).v)
% keyboard
% end

