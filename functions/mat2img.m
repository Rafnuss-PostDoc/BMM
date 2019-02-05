

folder ='Density_estimationMap_ImageOverlay/';
load('data/Density_estimationMap','g')
data = log(g.dens_est);
min_d=-2;
max_d=5;

% folder ='Density_simulationMap_ImageOverlay/';
% load('data/Density_simulationMap_reassemble','real_dens')
% data = log(real_dens);
% min_d=-5;
% max_d=5;

% folder ='FlightSpeed_estimationMap_ImageOverlay/';
% load('data/FlightSpeed_estimationMap','g')
% data = g.dens_est;

% folder ='FlightSpeed_simulationMap_ImageOverlay/';
% load('data/FlightSpeed_simulationMap_reassemble')
% data = real_dens;

% folder ='FlightDir_estimationMap_ImageOverlay/';
% load('data/FlightDir_estimationMap','gdir')
% data = mod(gdir.dd_est,360);

% folder ='FlightDir_simulationMap_ImageOverlay/';
% load('data/FlightDir_simulationMap_reassemble')
% data = mod(real_dir,360);

% folder = 'SinkSource_estimationMap_ImageOverlay/';
% data=W;
% min_d=-10; max_d=10;

% folder = 'SinkSource_estimationMap_ImageOverlay/';
% data=W;
% min_d=-10; max_d=10;

% folder = 'MTR_estimationMap_ImageOverlay/';
load('data/Density_estimationMap','g')
c = (g.dens_est);
load('data/FlightSpeed_estimationMap','g')
speed = g.dens_est;
% data = c .* speed * 60*60/1000;
% min_d=0; max_d=3000;

% histogram(data(:));

data = (data-min_d)/ (max_d-min_d);
data(data<0)=0;
data(data>1)=1;

cm=viridis(255);

A=ones(g.nlat,g.nlon);
% A=ones(g.nlat-2,g.nlon-2);

% Rain
% load total column rain water (kg/m^2)
% load('data/Rain_grid.mat');
% thr=.005;
% mean(TCRW(:)>thr);
% folder ='rain/';
% rgb=imresize(reshape(repmat([77, 77, 255]/255,g.nl,1),[g.nlat, g.nlon,3]),3,'nearest');

for i=1:size(data,3)
    
    img = flipud(data(:,:,i));

    % % Rain
    % img_rain = flipud(TCRW(:,:,i))>thr;

    img_nan=isnan(img);
    img(img_nan)=0;
    % rgb(img_nan(:),:)=0.8;

    rgb = imresize(reshape(cm(uint8(img(:)*255+1),:),[size(img),3]),3,'nearest');
    
    A2=A;
    % A2(img_nan)=0.8;
    % A2(img_rain)=1;
    A2(img_nan)=0;
   
    % imshow(rgb)
    filename = ['./figure/' folder datestr(g.time(i),'yyyy-mm-dd-HH-MM')];
    imwrite(rgb,[filename '.png'],'Alpha', imresize(A2,3,'nearest'))
    
    status = system(['gdal_translate -of Gtiff -a_ullr ' num2str([g.lon(1) g.lat(end) g.lon(end) g.lat(1)]) ' -a_srs EPSG:4326 ' filename '.png ' filename '_4326.tiff']);  assert(status==0)
    status = system(['gdalwarp -s_srs EPSG:4326 -t_srs EPSG:3857 -ts ' num2str([size(rgb,1) size(rgb,2)]) '  ' filename '_4326.tiff ' filename '_3857.tiff' ]); assert(status==0)
    status = system(['gdal_translate -of png ' filename '_3857.tiff ' filename '_3857.png']);  assert(status==0)
end
cd(['./figure/' folder])
status = system('del *.xml'); assert(status==0)
status = system('del *.tiff'); assert(status==0)
status = system('mkdir tmp'); assert(status==0)
status = system('move *3857.png ./tmp/'); assert(status==0)
status = system('del *.png'); assert(status==0)
cd(['./tmp/']); status = system('move *.png ..'); assert(status==0); cd ..
status = system('rmdir tmp'); assert(status==0)
cd ../..
