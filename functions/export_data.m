
load('data/Density_estimationMap' ,'g');

T=table();
T.date=g.time';


for i_lat=1:10:g.nlat
    for i_lon=1:10:g.nlon
        a=reshape(nanmean(nanmean(g.dens_sig,1),2),g.nt,[]);
        a=reshape(g.dens_sig(i_lat,i_lon,:,:),g.nt,[]);
        T.m3sigma=a(:,1);
        T.m2sigma=a(:,2);
        T.m1sigma=a(:,3);
        T.est=a(:,4);
        T.p1sigma=a(:,5);
        T.p2sigma=a(:,6);
        T.p3sigma=a(:,7);
        csvwrite(['figure/Density_latlon_csv/' num2str(g.lat(i_lat)) '_' num2str(g.lon(i_lon)) '.csv'], )
    end
end
writetable(T,['figure/Density_latlon_csv/all.csv'])


load('data/FlightSpeed_estimationMap','g')
data_est(:,3) = g.dens_est;
data_est(:,4) = g.dens_sig;
load('data/FlightDir_estimationMap','gdir')
data_est(:,5) = mod(gdir.dd_est,360);
data_est(:,6) = gdir.dd_sig;


load('data/Density_simulationMap','real_dens')
data_sim=(g.nlat*g.nlon*g.nat,3);
data_sim(:,1) = real_dens;
load('data/FlightSpeed_simulationMap_reassemble','real_dens')
data_sim(:,2) = real_dens;
load('data/FlightDir_simulationMap_reassemble','real_dir')
data_sim(:,3) = mod(real_dir,360);