
function quiver2img(u,v,rzd,g,folder)


u_isnan=single(~isnan(u));
v_isnan=single(~isnan(v));

u(isnan(u))=0;
v(isnan(v))=0;

lat2D_res=imresize(g.lat2D,rzd);
lon2D_res=imresize(g.lon2D,rzd);

figure('position',[0 0 900 1000])
for i=1627:size(u,3)
    if sum(sum(u(:,:,i)))==0
        copyfile('./BMM_web/blank.png',['./BMM_web/' folder datestr(g.time(i),'yyyy-mm-dd-HH-MM') '_3857.png'])
    else
        u_res = imresize(u(:,:,i),rzd);
        u_isnan_res = imresize(u_isnan(:,:,i),rzd);
        u_res=u_res./u_isnan_res;
        u_res(u_isnan_res<0.5)=nan;

        v_res = imresize(v(:,:,i),rzd);
        v_isnan_res = imresize(v_isnan(:,:,i),rzd);
        v_res=v_res./v_isnan_res;
        v_res(v_isnan_res<0.5)=nan;

        quiver(lon2D_res,lat2D_res,u_res/10,v_res/10,'AutoScale','off','linewidth',0.7,'color','k')
        % set(gca,'Visible','off')
        set(gca,'ytick',[]); set(gca,'xtick',[])
        set(gca, 'Color', 'none');
        axis([g.lon(1) g.lon(end) g.lat(1) g.lat(end)])
        filename = ['./BMM_web/' folder datestr(g.time(i),'yyyy-mm-dd-HH-MM')];
        export_fig([filename '.png'],'-transparent','-native','-m3')


        status = system(['gdal_translate -of Gtiff -a_ullr ' num2str([g.lon(1) g.lat(end) g.lon(end) g.lat(1)]) ' -a_srs EPSG:4326 ' filename '.png ' filename '_4326.tiff']);  assert(status==0)
        status = system(['gdalwarp -s_srs EPSG:4326 -t_srs EPSG:3857 -ts ' num2str(30*[size(v_res,1) size(v_res,2)]) '  ' filename '_4326.tiff ' filename '_3857.tiff' ]); assert(status==0)
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

end