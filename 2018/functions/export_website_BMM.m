clear all;
load('data/Density_estimationMap','g','gd');
load('data/Flight_estimationMap','guv');


dl=diff(g.lat(1:2))/2;

FC={};
FC.type="FeatureCollection";
FC.features=[];
F={};
F.type='Feature';
F.geometry={};
F.properties={};
F.geometry.type="Polygon";

tmplat = bsxfun(@plus, g.lat2D(g.latlonmask), [dl -dl -dl dl dl]);
tmplon = bsxfun(@plus, g.lon2D(g.latlonmask), [-dl -dl dl dl -dl]);

for i_l = 1:g.nlm
    F.geometry.coordinates(1,:,:)= [tmplon(i_l,:)' tmplat(i_l,:)'];
    F.properties.value=0;
    FC.features=[FC.features;F];
end

fileID = fopen('data/grid.geojson','w');
fwrite(fileID,jsonencode(FC));
fclose(fileID)


%% Export time for wesbite
time = datetime('01-Jan-2018 00:00'):1/24/4:datetime('01-Jan-2019 07:30');
[~,Locb] = ismember(time,g.time);

fileID = fopen('data/date.json','w');
fwrite(fileID,jsonencode(Locb));
fclose(fileID)

%% Export quiver point map

FC={};
FC.type="FeatureCollection";
FC.features=[];
F={};
F.type='Feature';
F.geometry={};
F.properties={};
F.geometry.type="Point";

tmp=[g.lon2D(g.latlonmask) g.lat2D(g.latlonmask)];
zoom=nan(g.nlat, g.nlon);
ii=[1 2 4 8 16 32 64 128];
for i=1:numel(ii)
    zoom(1:ii(i):end,1:ii(i):end)=numel(ii)-i;
end
zoom(~g.latlonmask)=nan;
imagesc(zoom<7); colorbar
zoom = zoom(g.latlonmask);

for i=1:size(tmp,1)
    F.properties.angle=round(rand(1)*365);
    F.properties.size=round(rand(1),2);
    F.properties.zoom=zoom(i);
    F.geometry.coordinates = tmp(i,:);
    FC.features=[FC.features;F];
end

fileID = fopen('data/quiver.geojson','w');
fwrite(fileID,jsonencode(FC));
fclose(fileID)

%% Export dens+flight bin

gd_dens_est = round(gd.dens_est*100);
guv_u_est = round((guv.u_est+30)*100);
guv_v_est = round((guv.v_est+30)*100);

% for i_ll=1:g.nlm
%     fileID = fopen(['data/bin/ll_' num2str(i_ll) '.bin'],'w');
%     fwrite(fileID,[gd_dens_est(i_ll,:)' guv_u_est(i_ll,:)' guv_v_est(i_ll,:)'],'uint16');
%     fclose(fileID);
% end
% 
% for i_t=1:g.nt
%     fileID = fopen(['data/bin/density_' num2str(i_t) '.bin'],'w');
%     fwrite(fileID,[gd_dens_est(:,i_t) guv_u_est(:,i_t) guv_v_est(:,i_t)],'uint16');
%     fclose(fileID);
% end

A = [gd_dens_est guv_u_est guv_v_est];
V1 = diff(int16(A(:)));
V = bitxor(bitshift(V1,1),-bitshift(V1,-15));
fileID = fopen('data/density.bin','w');
fwrite(fileID,[0; V1],'int16');
fclose(fileID);
! wsl gzip data/density.bin -f


% fwrite(fileID,[0; V],'int16');
%! wsl brotli data/bin/density.bin -f

%% Export radar data
load('data/Density_inference.mat','data');
FC={};
FC.type="FeatureCollection";
FC.features=[];
F={};
F.type='Feature';
F.geometry={};
F.properties={};
F.geometry.type="Point";

id = ismember(data.time,g.time);
id2 = ismember(g.time,data.time(id));

for i=1:data.nrad
    F.geometry.coordinates = [data.lon(i) data.lat(i)];
    F.properties.dens=zeros(g.nt,1);
    F.properties.dens(id2) =  round( data.dens_m(id,i)*100 );
    F.properties.dens(isnan(F.properties.dens))=0;
    F.properties.u=zeros(g.nt,1);
    F.properties.u(id2) = round( (data.us(id,i)+30)*100 );
    F.properties.u(isnan(F.properties.u))=0;
    F.properties.v=zeros(g.nt,1);
    F.properties.v(id2) = round( (data.vs(id,i)+30)*100 );
    F.properties.v(isnan(F.properties.v))=0;
    F.properties.name = data.name{i};
    FC.features=[FC.features;F];
end

fileID = fopen('data/radar.geojson','w');
fwrite(fileID,jsonencode(FC));
fclose(fileID);
 
%% Export rain

g_rain = reshape(g.rain(repmat(g.latlonmask,1,1,g.nt))>data.mask_rain_thr, g.nlm,g.nt);

fileID = fopen('data/rain.bin','w');
fwrite(fileID,g_rain);
fclose(fileID);
! wsl gzip data/rain.bin -f


%% Export Total MTR bird/km^2 * 
% m/s -> km/hr
dy = lldistkm([g.lat(1) g.lon(1)],[g.lat(2) g.lon(1)]);
dx = lldistkm([g.lat2D(:,1) g.lon2D(:,1)],[g.lat2D(:,1) g.lon2D(:,2)]);
area = repmat(dx*dy,1,g.nlon);
dt = 15/60;
MTR = gd.dens_est .* repmat(area(g.latlonmask),1,g.nt) .* sqrt(guv.u_est.^2 + guv.v_est.^2)/1000*60*60 .* dt;
TS_MTR = nansum(MTR,1);

fileID = fopen('data/MTR.json','w');
fwrite(fileID,jsonencode(round(TS_MTR/1000000)));
fclose(fileID)


