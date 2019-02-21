datestr(g.time,'yyyy-mm-dd-HH-MM')
load('data/Density_estimationMap','g')

load('data/Flight_estimationMap','guv')

load('data/Rain_grid.mat','TCRW');

geoidR = makerefmat('RasterSize', [g.nlat g.nlon], 'Latlim', [g.lat(1) g.lat(end)], 'Lonlim', [g.lon(1) g.lon(end)]);

for i=1:g.nt

    data = cat(3,g.dens_est(:,:,i),g.dens_q10(:,:,i),g.dens_q90(:,:,i),guv.u_est(:,:,i),guv.v_est(:,:,i),TCRW(:,:,i));
   
    rgb=data(:,:,1);
    filename = [folder datestr(g.time(i),'yyyy-mm-dd-HH-MM')];
    geotiffwrite([filename '.tif'], single(data), geoidR)
end
cd([folder])
status = system('del *.png'); assert(status==0)
status = system('del *_4326.tiff'); assert(status==0)
cd ../..
