% function SunriseSunset(d, start_date,end_date)

% Data 
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



% Grid 
lat_d=43:68; lat_d(lat_d==0)=0.0000001;
lon_d=-5:30;lon_d(lon_d==0)=0.0000001;
tim_d = start_date-1/2:end_date+1/2;
[LAT,LON,TIM]=ndgrid(lat_d,lon_d,tim_d);

sunrs_e=nan(numel(lat_d),numel(lon_d),numel(tim_d));
sunrs_b=nan(numel(lat_d),numel(lon_d),numel(tim_d));
for i_lat=1:5:numel(lat_d)
    for i_lon=1:5:numel(lon_d)
        for i_t=1:5:numel(tim_d)
%             dd = webread(['https://api.sunrise-sunset.org/json?lat=' num2str(lat_d(i_lat)) '&lng=' num2str(lon_d(i_lon)) '&date=' datestr(tim_d(i_t),'yyyy-mm-dd')]);
%             assert(strcmp(dd.status,'OK'))
            sunrs_b(i_lat,i_lon,i_t) = a(i_lat,i_lon,i_t);% datenum(string(dd.results.civil_twilight_begin));
            sunrs_e(i_lat,i_lon,i_t) = b(i_lat,i_lon,i_t);% datenum(string(dd.results.civil_twilight_end));
            % pause(5)
        end
    end
end

in = isnan(sunrs_e);

F=griddedInterpolant(reshape(LAT(~in),6,8,5),reshape(LON(~in),6,8,5),datenum(reshape(TIM(~in),6,8,5)),reshape(sunrs_e(~in),6,8,5));
sunrs_e=F(LAT,LON,datenum(TIM));

F=griddedInterpolant(reshape(LAT(~in),6,8,5),reshape(LON(~in),6,8,5),datenum(reshape(TIM(~in),6,8,5)),reshape(sunrs_b(~in),6,8,5));
sunrs_b=F(LAT,LON,datenum(TIM));

sunrs_b = repmat(reshape(datenum(tim_d-1/2),1,1,[]),size(sunrs_b,1),size(sunrs_b,2),1)+sunrs_b-datenum('01-01-2018');
sunrs_e = repmat(reshape(datenum(tim_d-1/2),1,1,[]),size(sunrs_e,1),size(sunrs_e,2),1)+sunrs_e-datenum('01-01-2018');

save('./data/sunrisesunset_grid.mat','lat_d','lon_d','tim_d','sunrs_e','sunrs_b')
% load('./data/sunrisesunset_grid.mat')












% OLD METHOD, not computing civil twilight
% 
% Sun elevation
% for i_d=1:numel(d)
%     tmp1 = struct('latitude',d(i_d).lat,'longitude',d(i_d).lon,'altitude',d(i_d).height);
%     tmp2 = struct('UTCOffset', 1,'year', year(d(i_d).time), 'month', month(d(i_d).time),'day', day(d(i_d).time),'hour', hour(d(i_d).time), 'minute', minute(d(i_d).time),'second', second(d(i_d).time));
%     [~, d(i_d).sun, ~]=pvl_spa_perso(tmp2, tmp1);
% end

%Sunset and rise
% for i_d=1:numel(d)
%     [tmp1, tmp2] = SunriseSunset(d(i_d).lon, d(i_d).lat, 1);
%     tmp3 = datenum(datestr(d(i_d).time,'yyyy-mm-dd'))-datenum(2016,1,0,0,0,0);
% 
%     d(i_d).sunrise=tmp1(tmp3)';
%     d(i_d).sunset=tmp2(tmp3)';
%     
%     d(i_d).sunrise(hour(d(i_d).time)>12) = tmp1(tmp3(hour(d(i_d).time)>12)+1);
%     d(i_d).sunset(hour(d(i_d).time)<12) = tmp2(tmp3(hour(d(i_d).time)<12)-1);
%   
%     d(i_d).sunset_time = datetime([datestr(d(i_d).time,'dd-mmm-yyyy') datestr(d(i_d).sunset/24)]);
%     d(i_d).sunrise_time = datetime([datestr(d(i_d).time,'dd-mmm-yyyy') datestr(d(i_d).sunrise/24)]);
% end
% 
% function [sunrise,sunset] = SunriseSunset(long, lat, UTCoff)
% longCorr = 4*(long - 15*UTCoff);                         % longitudinal correction
% days = 1:365;
% B = 360*(days - 81)/365;
% EoTCorr = 9.87*sind(2*B) - 7.53*cosd(B) - 1.5*sind(B);  % Equation of Time correction
% solarCorr = longCorr + EoTCorr;
% 
% delta = asind(sind(23.45)*sind(360*(days - 81)/365));    % Solar declination
% 
% sunrise = 12 - acosd(-tand(lat)*tand(delta))/15 - solarCorr/60;
% sunset  = 12 + acosd(-tand(lat)*tand(delta))/15 - solarCorr/60;
% end