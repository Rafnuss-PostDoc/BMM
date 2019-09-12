
addpath('./functions/')

%% Export raw data (Zenodo)
load('../1-Cleaning/data/dc_corr.mat'); 
for i_d=1:numel(dc)
    dc(i_d).alt=dc(i_d).interval*(1:dc(i_d).levels)-dc(i_d).interval/2;
end

dc_zenodo = rmfield(dc,[{'date','stime','etime','n_all','DBZH','w','ff','dd','NNT','interval','levels', 'cc', 't', 'crwc'} q.s]); 

for i_d=1:numel(dc_zenodo)
    
    dc_zenodo(i_d).u = round(dc_zenodo(i_d).u,2);
    dc_zenodo(i_d).v = round(dc_zenodo(i_d).v,2);
    dc_zenodo(i_d).dens = round(dc_zenodo(i_d).dens,2);
    dc_zenodo(i_d).us = round(dc_zenodo(i_d).us,2);
    dc_zenodo(i_d).vs = round(dc_zenodo(i_d).vs,2);
    dc_zenodo(i_d).denss = round(dc_zenodo(i_d).denss,2);
    
    fileID = fopen(['./zenodo/dc_' dc_zenodo(i_d).name '.json'],'w');
    fprintf(fileID,jsonencode(dc_zenodo(i_d)));
    fclose(fileID);
end

%% Export Interpolation Estimation with uncertainty (Zenodo)
load('../4-Estimation/data/Density_estimationMap');
load('../6-Flight/data/Flight_estimationMap');

id = ~isnan(gd.dens_est) & ~isnan(guv.u_est);
t=table();
t.density_estimation = round(gd.dens_est(id),2);
t.density_quantile10 = round(gd.dens_q10(id),2);
t.density_quantile90 = round(gd.dens_q90(id),2);
t.speedu_estimation = round(guv.u_est(id),2);
t.speedu_std = round(guv.u_sig(id),2);
t.speedv_estimation = round(guv.v_est(id),2);
t.speedv_std = round(guv.v_sig(id),2);
t.latitude = g.lat3D(id);
t.longitude = g.lon3D(id);
t.time = g.time3D(id);

writetable(t,'./zenodo/interpolation_result.csv')


%% Export for MangoDB
load('../4-Estimation/data/Density_estimationMap');
load('../6-Flight/data/Flight_estimationMap');

dlat = lldistkm([g.lat(1) g.lon(1)],[g.lat(2) g.lon(1)]);
dlon = lldistkm([g.lat2D(:,1) g.lon2D(:,1)],[g.lat2D(:,1) g.lon2D(:,2)]);
area = repmat(dlon*dlat,1,g.nlon);
area( ~g.mask_water | ~g.mask_distrad)=nan;
area = repmat(area,1,1,g.nt);
area(~g.mask_rain)=nan;

g_dens_est=gd.dens_est;
g_dens_est(~g.mask_rain)=nan;
g_dens_q10=gd.dens_q10;
g_dens_q10(~g.mask_rain)=nan;
g_dens_q90=gd.dens_q90;
g_dens_q90(~g.mask_rain)=nan;

a = reshape(nansum(nansum(area,1),2),g.nt,[]);
est = reshape( nansum( nansum( g_dens_est .* area ,1) ,2),[],1) ./ a;
idt=est~=0;
idt(idt==0 & ([idt(2:end);1]==1 | [1;idt(1:end-1)]==1))=1; % Adding a buffer of 1.
idt = find(idt);
        
dd=2;
i=1;
density=nan(numel(idt),3,g.nlm);
uv=nan(numel(idt),2,g.nlm);

latlon=nan(3,g.nlm);
for i_lat=1:dd:g.nlat-dd
    for i_lon=1:dd:g.nlon-dd
        if any(g.latlonmask(i_lat+(0:dd-1),i_lon+(0:dd-1)))
            surface=max(reshape(nansum(nansum(area(i_lat+(0:dd-1),i_lon+(0:dd-1),idt),1),2),[],1));
            latlon(:,i) = round([mean(g.lat(i_lat+(0:dd-1))) mean(g.lon(i_lon+(0:dd-1))) surface],2);
            
            density(:,1,i) = reshape( nanmean(nanmean( g_dens_est(i_lat+(0:dd-1),i_lon+(0:dd-1),idt),1),2),[],1);
            density(:,2,i) = reshape( nanmean(nanmean( g_dens_q10(i_lat+(0:dd-1),i_lon+(0:dd-1),idt),1),2),[],1);
            density(:,3,i) = reshape( nanmean(nanmean( g_dens_q90(i_lat+(0:dd-1),i_lon+(0:dd-1),idt),1),2),[],1);
            
            uv(:,1,i) = reshape( nanmean(nanmean( guv.u_est(i_lat+(0:dd-1),i_lon+(0:dd-1),idt),1),2),[],1);
            uv(:,2,i) = reshape( nanmean(nanmean( guv.v_est(i_lat+(0:dd-1),i_lon+(0:dd-1),idt),1),2),[],1);

            i=i+1;
        end
    end
end
density=density(:,:,1:i-1);
density(isnan(density))=0;
uv=uv(:,:,1:i-1);
uv(isnan(uv))=0;
latlon=latlon(:,1:i-1);

% fileID = fopen('BMM_web/exportDensityGrid.json','w');
% fprintf(fileID,jsonencode({latlon',round(permute(density,[3,2,1]),2),round(all_Density',2)}));
% fclose(fileID);

fileID = fopen('BMM_web/API/exportEst_time.json','w');
fprintf(fileID,jsonencode({datestr(g.time(idt),'yyyy-mm-dd HH:MM')}));
fclose(fileID);

fileID = fopen('BMM_web/API/exportEst_grid.json','w');
fprintf(fileID,jsonencode(latlon'));
fclose(fileID);

% fileID = fopen('BMM_web/exportEst_density_all.json','w');
% fprintf(fileID,jsonencode({round(all_Density',2)}));
% fclose(fileID);
% fileID = fopen('BMM_web/exportEst_density.json','w');
% fprintf(fileID,jsonencode({round(round(permute(density,[3,2,1]),2),2)}));
% fclose(fileID);
% fileID = fopen('BMM_web/exportEst_UV.json','w');
% fprintf(fileID,jsonencode({round(permute(uv,[3,1,2]),2)}));
% fclose(fileID);

% Mangodb
clear d
for i=1:size(density,3)
    d(i).density.est = round(density(:,1,i),2);
    d(i).density.q10 = round(density(:,2,i),2);
    d(i).density.q90 = round(density(:,3,i),2);
    d(i).u = round(uv(:,1,i),2);
    d(i).v = round(uv(:,2,i),2);
    d(i).id = i;
end
str = jsonencode(d);
str=strrep(str(2:end-1),'},{','}\n{');
str=strrep(str,'id','_id');
fileID = fopen('BMM_web/API/exportEst_mongodb.json','w');
fprintf(fileID,str);
fclose(fileID);


%% Simulation
load('../5-Simulation/data/Density_simulationMap_reassemble_ll');

% Compute area
dlat = lldistkm([g.lat(1) g.lon(1)],[g.lat(2) g.lon(1)]);
dlon = lldistkm([g.lat2D(:,1) g.lon2D(:,1)],[g.lat2D(:,1) g.lon2D(:,2)]);
area = repmat(dlon*dlat,1,g.nlon);
area( ~g.mask_water | ~g.mask_distrad)=nan;
area = repmat(area,1,1,g.nt);
area(~g.mask_rain)=nan;

% Date
a = reshape(nansum(nansum(area,1),2),g.nt,[]);
est = reshape( nansum( nansum( gd.dens_est .* area ,1) ,2),[],1) ./ a;
idt=est~=0;
idt2 = find(idt);
idt3 = idt2(1:2:end);
idtb = find([diff(idt) ;1] ~= 0 | [1; diff(idt)] ~= 0); % read buffer
idt = unique([idt3;idtb]);

% subsample real
%[~,real_sort]=sort(reshape(nansum(nansum(real_nb_ll),2),1,[]));
real_dens_ll_s = real_dens_ll(:,idt,1:50);
clear real_dens_ll

% Compute number of bird accounting for area and rain
real_nb_ll = real_dens_ll_s .* repmat(reshape(area(repmat(g.latlonmask,1,1,numel(idt))),g.nlm,[]),1,1,size(real_dens_ll_s,3));
clear real_dens_ll_s


% Determine latlon coordinate for masked
mask_latlon = [reshape(g.lat2D(g.latlonmask),[],1) reshape(g.lon2D(g.latlonmask),[],1)];

dd=4;
i=1;
sim=nan(numel(idt),size(real_nb_ll,3),g.nlm);

latlon=nan(3,g.nlm);
for i_lat=1:dd:g.nlat-dd
    for i_lon=1:dd:g.nlon-dd
        if any(g.latlonmask(i_lat+(0:dd-1),i_lon+(0:dd-1)))
            surface=max(reshape(nansum(nansum(area(i_lat+(0:dd-1),i_lon+(0:dd-1),idt),1),2),[],1));
            latlon(:,i) = round([mean(g.lat(i_lat+(0:dd-1))) mean(g.lon(i_lon+(0:dd-1))) surface],2);

            id=ismember(mask_latlon,[reshape(g.lat2D(i_lat+(0:dd-1),i_lon+(0:dd-1)),[],1) reshape(g.lon2D(i_lat+(0:dd-1),i_lon+(0:dd-1)),[],1)],'rows');
            
            sim(:,:,i) = permute(nansum(real_nb_ll(id,:,:)),[2 3 1]);

            i=i+1;
        end
    end
end
sim=sim(:,:,1:i-1);
sim(isnan(sim))=0;
latlon=latlon(:,1:i-1);

% Writing files
fileID = fopen('BMM_web/API/exportSim_grid.json','w');
fprintf(fileID,jsonencode(latlon'));
fclose(fileID);

fileID = fopen('BMM_web/API/exportSim_time.json','w');
fprintf(fileID,jsonencode({datestr(g.time(idt),'yyyy-mm-dd HH:MM')}));
fclose(fileID);


[I,J,K] = ndgrid(1:size(sim,1),1:size(sim,2),1:size(sim,3));

d = struct('sim', num2cell(sim(:)), 'loc_id', num2cell(K(:)), 'real_id', num2cell(J(:)), 'time_id', num2cell(I(:)));

str = jsonencode(d);
str=strrep(str(2:end-1),'},{','}\n{');
fileID = fopen('BMM_web/API/exportSim_mongodb.json','w');
fprintf(fileID,str);
fclose(fileID);


%% Binary

% fileID = fopen('BMM_web/exportDensityGrid.bin','w');
% fwrite(fileID,cat(3,density,all_Density),'single');
% fclose(fileID);
% 
% fileID = fopen('BMM_web/exportDensityGrid_latlon_time.json','w');
% fprintf(fileID,jsonencode({latlon, datestr(g.time(idt),'yyyy-mm-dd HH:MM')}));
% fclose(fileID);

%% Export ImageOverlay

% Estimation 
folder ='Density_estimationMap_ImageOverlay/';
data = log10(gd.dens_est);
min_d=0;
max_d=2;
cm=viridis(255); % https://ch.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps
mat2img(data,min_d,max_d,cm,g,folder,0)

% Simulation
folder ='Density_simulationMap_ImageOverlay/';
i_real=1;
tmp = nan(g.nlat,g.nlon,g.nt);
tmp(repmat(g.latlonmask,1,1,g.nt)) = real_dens_ll(:,:,i_real);
data = log10(tmp);
min_d=0;
max_d=2;
mat2img(data,min_d,max_d,cm,g,folder,0)

% Sink/source
% folder = 'SinkSource_estimationMap_ImageOverlay/';
% load('../8-SinkSource/data/SinkSourceSim.mat','W')
% data=W;
% min_d=-15; max_d=15;
% cm=colormap(brewermap([],'Spectral'));
% cm = interp1(linspace(1,255,64),cm,1:255);
% mat2img(data,min_d,max_d,g,folder)

% Rain
load('../1-Cleaning/data/Rain_grid.mat'); % load total column rain water (kg/m^2)
mean(TCRW(:)>g.mask_rain_thr);
folder ='rain2/';
data = TCRW;
% 
mat2img(data,min_d,max_d,cm,g,folder,1)


%% Export Quiver

load('../6-Flight/data/Flight_estimationMap')
folder='Quiver_est2/';
u=guv.u_est;
v=guv.v_est;
rzd =1/4;
quiver2img(u,v,rzd,g,folder)

load('../6-Flight/data/Flight_simulationMap_real_ll')
folder='Quiver_sim2/';
u = nan(g.nlat,g.nlon,g.nt);
v = nan(g.nlat,g.nlon,g.nt);
u(repmat(g.latlonmask,1,1,g.nt)) = real_u_ll(:,:,1);
v(repmat(g.latlonmask,1,1,g.nt)) = real_v_ll(:,:,1);
rzd =1/4;
quiver2img(u,v,rzd,g,folder)

%% OLD NOT USED
% % Export old csv
% T=table();
% T.date=g.time';
% for i_lat=1:10:g.nlat
%     for i_lon=1:10:g.nlon
%         a=reshape(nanmean(nanmean(g.dens_sig,1),2),g.nt,[]);
%         a=reshape(g.dens_sig(i_lat,i_lon,:,:),g.nt,[]);
%         T.m3sigma=a(:,1);
%         T.m2sigma=a(:,2);
%         T.m1sigma=a(:,3);
%         T.est=a(:,4);
%         T.p1sigma=a(:,5);
%         T.p2sigma=a(:,6);
%         T.p3sigma=a(:,7);
%         csvwrite(['figure/Density_latlon_csv/' num2str(g.lat(i_lat)) '_' num2str(g.lon(i_lon)) '.csv'], )
%     end
% end
% writetable(T,['figure/Density_latlon_csv/all.csv'])


% load('data/FlightSpeed_estimationMap','g')
% data_est(:,3) = g.dens_est;
% data_est(:,4) = g.dens_sig;
% load('data/FlightDir_estimationMap','gdir')
% data_est(:,5) = mod(gdir.dd_est,360);
% data_est(:,6) = gdir.dd_sig;
% 
% 
% load('data/Density_simulationMap','real_dens')
% data_sim=(g.nlat*g.nlon*g.nat,3);
% data_sim(:,1) = real_dens;
% load('data/FlightSpeed_simulationMap_reassemble','real_dens')
% data_sim(:,2) = real_dens;
% load('data/FlightDir_simulationMap_reassemble','real_dir')
% data_sim(:,3) = mod(real_dir,360);