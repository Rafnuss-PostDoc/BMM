
% folder ='Density_estimationMap_ImageOverlay/';
% load('data/Density_estimationMap')
% data = log10(gd.dens_est);
% min_d=0;
% max_d=2;

folder = 'SinkSource_estimationMap_ImageOverlay/';
load('data/SinkSourceSim.mat','W')
data=W;
min_d=-15; max_d=15;
gold=g;
g = gext;
g.time=gold.time;

% folder ='Density_simulationMap_ImageOverlay/';
% load('data/Density_simulationMap')
% data = log10(real_dens(:,:,:,i_real));
% min_d=0;
% max_d=3;

% histogram(data(:));

data = (data-min_d)/ (max_d-min_d);
data(data<0)=0;
data(data>1)=1;

cm=viridis(255);
% cm=colormap(brewermap([],'Spectral'))
% cm = interp1(linspace(1,255,64),cm,1:255);

% Rain
% load total column rain water (kg/m^2)
% load('data/Rain_grid.mat');
% mean(TCRW(:)>g.mask_rain_thr);
% folder ='rain/';
% data = TCRW;
% rgb=imresize(reshape(repmat([77, 77, 255]/255,g.nl,1),[g.nlat, g.nlon,3]),3,'nearest');


for i=1:size(data,3)
    if (sum(sum(data(:,:,i)))==0 || all(all(isnan(data(:,:,i)))) )
        copyfile('./BMM_web/blank_sinksource.png',['./BMM_web/' folder datestr(g.time(i),'yyyy-mm-dd-HH-MM') '_3857.png'])
    else
        img = flipud(data(:,:,i));
        img_nan=isnan(img);
        img(img_nan)=0;
        rgb = imresize(reshape(cm(uint8(img(:)*255+1),:),[size(img),3]),3,'nearest');
        A=ones(g.nlat,g.nlon);
        A(img_nan)=0;
        
        % Rain
%         img_rain = flipud(data(:,:,i))>g.mask_rain_thr;
%         A=zeros(g.nlat,g.nlon);
%         A(img_rain)=1;

        % imshow(rgb); drawnow;
        filename = ['./BMM_web/' folder datestr(g.time(i),'yyyy-mm-dd-HH-MM')];
        imwrite(rgb,[filename '.png'],'Alpha', imresize(A,3,'nearest'))

        status = system(['gdal_translate -of Gtiff -a_ullr ' num2str([g.lon(1) g.lat(end) g.lon(end) g.lat(1)]) ' -a_srs EPSG:4326 ' filename '.png ' filename '_4326.tiff']);  assert(status==0)
        status = system(['gdalwarp -s_srs EPSG:4326 -t_srs EPSG:3857 -ts ' num2str([size(rgb,1) size(rgb,2)]) '  ' filename '_4326.tiff ' filename '_3857.tiff' ]); assert(status==0)
        status = system(['gdal_translate -of png ' filename '_3857.tiff ' filename '_3857.png']);  assert(status==0)
    end
end

cd(['./BMM_web/' folder])
status = system('del *.xml'); assert(status==0)
status = system('del *.tiff'); assert(status==0)
status = system('mkdir tmp'); assert(status==0)
status = system('move *3857.png ./tmp/'); assert(status==0)
status = system('del *.png'); assert(status==0)
cd(['./tmp/']); status = system('move *.png ..'); assert(status==0); cd ..
status = system('rmdir tmp'); assert(status==0)
cd ../..