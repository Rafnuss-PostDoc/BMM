
load('data/FlightSpeed_estimationMap','g')
speed = g.dens_est;
load('data/FlightDir_estimationMap','gdir')
dir = mod(gdir.dd_est,360);
folder='Quiver_est_4/';

load('data/FlightSpeed_simulationMap_reassemble','real_dens')
speed = real_dens;
load('data/FlightDir_simulationMap_reassemble','real_dir')
dir = mod(real_dir,360);
folder='Quiver_sim_4/';

dir=mod(90-dir,360);
deltalat=speed.*cosd(dir);
deltalon=speed.*sind(dir);

% worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
% setm(gca,'Frame','off');
% setm(gca,'grid','off')
% load coastlines;
% plotm(coastlat, coastlon,'k')

rzd=1/4;

deltalat_a=deltalat;
deltalat_a(isnan(deltalat_a))=0;
deltalon_a=deltalon;
deltalon_a(isnan(deltalon_a))=0;

deltalat_b=single(~isnan(deltalat));
deltalon_b=single(~isnan(deltalon));

lat2D_res=imresize(g.lat2D,rzd);
lon2D_res=imresize(g.lon2D,rzd);


figure('position',[0 0 900 1000])

for i=1:size(deltalat_a,3)
  
    deltalat_a_res = imresize(deltalat_a(:,:,i),rzd);
    deltalat_b_res = imresize(deltalat_b(:,:,i),rzd);
    deltalat_res=deltalat_a_res./deltalat_b_res;
    deltalat_res(deltalat_b_res<0.5)=nan;
        
    deltalon_a_res = imresize(deltalon_a(:,:,i),rzd);
    deltalon_b_res = imresize(deltalon_b(:,:,i),rzd);
    deltalon_res=deltalon_a_res./deltalon_b_res;
    deltalon_res(deltalon_b_res<0.5)=nan;
    
    quiver(lon2D_res,lat2D_res,deltalat_res,deltalon_res,'linewidth',0.7,'color','k')
    % set(gca,'Visible','off')
    set(gca,'ytick',[]); set(gca,'xtick',[])
    set(gca, 'Color', 'none');
    axis([g.lon(1) g.lon(end) g.lat(1) g.lat(end)])
    filename = ['./figure/' folder datestr(g.time(i),'yyyy-mm-dd-HH-MM')];
    export_fig([filename '.png'],'-transparent','-native','-m3')
    
   
    status = system(['gdal_translate -of Gtiff -a_ullr ' num2str([g.lon(1) g.lat(end) g.lon(end) g.lat(1)]) ' -a_srs EPSG:4326 ' filename '.png ' filename '_4326.tiff']);  assert(status==0)
    status = system(['gdalwarp -s_srs EPSG:4326 -t_srs EPSG:3857 -ts ' num2str(30*[size(deltalon_res,1) size(deltalon_res,2)]) '  ' filename '_4326.tiff ' filename '_3857.tiff' ]); assert(status==0)
    status = system(['gdal_translate -of png ' filename '_3857.tiff ' filename '_3857.png']);  assert(status==0)
    
    % clear(h);
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