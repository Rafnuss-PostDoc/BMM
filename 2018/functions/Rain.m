function Rain(d, start_date,end_date,g)

% Data 
file = ['../ECMWF/output_20160919-20161010_srf_europe.nc'];
tcrw = ncread(file,'tcrw');
lat = double(ncread(file,'latitude'));
lon = double(ncread(file,'longitude'));
time = datetime(single(ncread(file,'time'))/24+datenum('1900-01-01'), 'ConvertFrom', 'datenum');

F = griddedInterpolant({-lat,lon,datenum(time)},permute(tcrw,[2 1 3]),'spline');

figure; hold on;
imagesc(-lat,lon,tcrw(:,:,30))
scatter(-g.lat2D(:),g.lon2D(:),'.k')


TCRW = F({-g.lat,g.lon,datenum(g.time)});

save('data/Rain_grid.mat','TCRW')

end
