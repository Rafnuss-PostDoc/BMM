script


load('data/dens_corr','dc');

%% Get elevatation around

files={'gt30w020n90'};%,'gt30e020n90'}; % ,'gt30e020n40','gt30w020n40'
[DEM,GeoCellRef] = geotiffread(['C:\Users\rnussba1\Documents\MATLAB\bin\' files{1}]);
GeoCellRef_lon = linspace(GeoCellRef.LongitudeLimits(1),GeoCellRef.LongitudeLimits(2),GeoCellRef.RasterSize(2));
GeoCellRef_lat = linspace(GeoCellRef.LatitudeLimits(2),GeoCellRef.LatitudeLimits(1),GeoCellRef.RasterSize(1));
[GeoCellRef_LON,GeoCellRef_LAT] = meshgrid(GeoCellRef_lon,GeoCellRef_lat);

figure('position',[0 0 800 600]);
worldmap([floor(min([dc.lat])) ceil(max([dc.lat]))], [floor(min([dc.lon])) ceil(max([dc.lon]))]);
geoshow(DEM,GeoCellRef,'DisplayType','texturemap'); demcmap(DEM)

% Set water level to 0
DEM(DEM<100)=0;

% projection angle
alpha=-90:10:90;

% DLAT, DLON
dlat = lldistkm([dc.lat; dc.lon]', [dc.lat; dc.lon]' + repmat([GeoCellRef.CellExtentInLatitude 0],numel(dc),1));
dlon = lldistkm([dc.lat; dc.lon]', [dc.lat; dc.lon]' + repmat([0 GeoCellRef.CellExtentInLongitude],numel(dc),1));

% volume in km^3
V = nan(numel(dc),numel(alpha));


for i_d=1:numel(dc)
    % pre-select coordinate within 1° lat and lon to reduce computational
    % time of lldistkm.
    id=false(GeoCellRef.RasterSize);
    id(abs(GeoCellRef_lat-dc(i_d).lat)<1 , abs(GeoCellRef_lon-dc(i_d).lon)<1) = true;
    distkm = lldistkm([dc(i_d).lat dc(i_d).lon],[GeoCellRef_LAT(id) GeoCellRef_LON(id)]);
    distkm_id = distkm<dc(i_d).maxrange;
    id_1 = find(id);
    id(id_1(~distkm_id))=0;
    
    % Find DEM of radar -> ELEVATION INCORRECT COMPARED TO GOOGLE MAP
    % ~0-20m (more for higher altitude)
%     [~,i_min] = min(distkm); 
%     d(i_d).heightDEM = DEM(id_1(i_min));

    % We imrotate only the matrix which contain value within distance.
    % we subsample to the rectangle encompasing the radar scan
    id_any_1 = any(id,2); id_any_2 = any(id);
    id_s = id(id_any_1, id_any_2);
    DEM_id = DEM(id_any_1, id_any_2);
    DEM_id(~id_s)=0;
    
    for i_alpha=1:numel(alpha)
        J = imrotate(id_s,alpha(i_alpha));
        DEM_id_r = imrotate(DEM_id,alpha(i_alpha));
        V(i_d,i_alpha) = sum( (dc(i_d).height-max(double(DEM_id_r)))/1000 .* sum(J)) *dlat(i_d)*dlon(i_d);
%         subplot(5,4,i_alpha); 
%         % imagesc(DEM_id_r); 
%         plot((d(i_d).height-max(double(DEM_id_r))))
%         title(['a=' num2str(alpha(i_alpha)) '   V=' num2str(V(i_d,i_alpha))])
    end
    dc(i_d).VolBelow=V(i_d,:);
%     
    figure; hold on; title(dc(i_d).name)
    %scatter3(GeoCellRef_LAT(:), GeoCellRef_LON(:), DEM(:),'.k')
    surf(GeoCellRef_LON(id_any_1, id_any_2), GeoCellRef_LAT(id_any_1, id_any_2), DEM_id,'EdgeColor','none')
    
    plot3(dc(i_d).lon, dc(i_d).lat, dc(i_d).height,'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10);
    xlabel('lat'); ylabel('lon'); zlabel('elev');
    view(3)
end

V2=V./repmat(pi*[dc.maxrange]'.^2*0.2,1,numel(alpha));

% figure;
figure; 
for i_d=1:numel(dc)
    subplot(6,7,i_d); 
    x1=deg2rad(alpha); y1=V2(i_d,:); x1(V2(i_d,:)<=0)=NaN; y1(V2(i_d,:)<=0)=NaN;
    polarplot(x1,abs(y1))
    hold on;
    x1=deg2rad(alpha); y1=V2(i_d,:); x1(V2(i_d,:)>=0)=NaN; y1(V2(i_d,:)>=0)=NaN;
    polarplot(x1,abs(y1));
    ax=gca;
    ax.ThetaDir = 'clockwise';
    ax.ThetaLim=[-90 90];
    ax.ThetaZeroLocation='top';
    rlim([0 1.2])
    %thetaticks('');
    title(dc(i_d).name)
end

% Find elevation with API OPEN_ELEVATION. NOT ALWAYS CORRECT. EXECEL
% spreadsheet has corrected value with https://www.freemaptools.com/elevation-finder.htm
strloc = join(strcat( strcat(cellstr(num2str([dc.lat]')),',') , cellstr(num2str([dc.lon]')) )','|');
strloc = replace(strloc{1},' ','');
tmp = webread(['https://api.open-elevation.com/api/v1/lookup?locations=' strloc ]); % num2str(d(i_d).lat) ',' num2str(d(i_d).lon)
tmp = jsondecode('{"results": [{"latitude": 51.1917, "elevation": -1, "longitude": 3.0642}, {"latitude": 49.9143, "elevation": 546, "longitude": 5.5056}, {"latitude": 49.6583, "elevation": 861, "longitude": 13.8178}, {"latitude": 49.5011, "elevation": 744, "longitude": 16.7885}, {"latitude": 53.564, "elevation": 1, "longitude": 6.7483}, {"latitude": 54.0044, "elevation": 93, "longitude": 10.0468}, {"latitude": 51.1246, "elevation": 225, "longitude": 13.7686}, {"latitude": 49.5407, "elevation": 774, "longitude": 12.4028}, {"latitude": 53.3394, "elevation": 4, "longitude": 7.025}, {"latitude": 51.4055, "elevation": 153, "longitude": 6.9669}, {"latitude": 51.3112, "elevation": 560, "longitude": 8.802}, {"latitude": 47.8736, "elevation": 1487, "longitude": 8.0036}, {"latitude": 52.4601, "elevation": 58, "longitude": 9.6945}, {"latitude": 48.0431, "elevation": 668, "longitude": 10.2204}, {"latitude": 50.5001, "elevation": 844, "longitude": 11.135}, {"latitude": 50.1097, "elevation": 553, "longitude": 6.5485}, {"latitude": 49.9859, "elevation": 204, "longitude": 8.714}, {"latitude": 52.6486, "elevation": 145, "longitude": 13.8578}, {"latitude": 54.1757, "elevation": 2, "longitude": 12.0581}, {"latitude": 48.1747, "elevation": 633, "longitude": 12.1018}, {"latitude": 48.5853, "elevation": 733, "longitude": 9.7828}, {"latitude": 52.1601, "elevation": 153, "longitude": 11.1761}, {"latitude": 55.1127, "elevation": 162, "longitude": 14.8875}, {"latitude": 55.1731, "elevation": 2, "longitude": 8.552}, {"latitude": 57.4893, "elevation": 85, "longitude": 10.1365}, {"latitude": 55.3262, "elevation": 31, "longitude": 12.4493}, {"latitude": 56.024, "elevation": 120, "longitude": 10.0246}, {"latitude": 59.3977, "elevation": 33, "longitude": 24.6021}, {"latitude": 58.4823, "elevation": 131, "longitude": 25.5187}, {"latitude": 36.8325, "elevation": 479, "longitude": -2.08222}, {"latitude": 39.4289, "elevation": 605, "longitude": -6.28528}, {"latitude": 41.4081, "elevation": 613, "longitude": 1.88472}, {"latitude": 43.1689, "elevation": 592, "longitude": -8.52694}, {"latitude": 41.9956, "elevation": 873, "longitude": -4.60278}, {"latitude": 28.0186, "elevation": 1760, "longitude": -15.6144}, {"latitude": 40.1758, "elevation": 702, "longitude": -3.71361}, {"latitude": 36.6133, "elevation": 1131, "longitude": -4.65917}, {"latitude": 38.2644, "elevation": 1187, "longitude": -1.18972}, {"latitude": 39.3797, "elevation": 109, "longitude": 2.785}, {"latitude": 43.4625, "elevation": 867, "longitude": -6.30194}, {"latitude": 37.6875, "elevation": 501, "longitude": -6.33444}, {"latitude": 43.4033, "elevation": 554, "longitude": -2.84194}, {"latitude": 39.1761, "elevation": 166, "longitude": -0.25211}, {"latitude": 41.7339, "elevation": 811, "longitude": -0.54583}, {"latitude": 60.9039, "elevation": 181, "longitude": 27.1081}, {"latitude": 61.7673, "elevation": 157, "longitude": 23.0764}, {"latitude": 61.907, "elevation": 186, "longitude": 29.7977}, {"latitude": 60.1285, "elevation": 113, "longitude": 21.6434}, {"latitude": 62.8626, "elevation": 107, "longitude": 27.3815}, {"latitude": 67.1391, "elevation": 0, "longitude": 26.8969}, {"latitude": 64.7749, "elevation": 93, "longitude": 26.3189}, {"latitude": 60.2706, "elevation": 160, "longitude": 24.869}, {"latitude": 63.1048, "elevation": 0, "longitude": 23.8209}, {"latitude": 50.1358, "elevation": 69, "longitude": 1.83472}, {"latitude": 42.1297, "elevation": 40, "longitude": 9.49639}, {"latitude": 50.1283, "elevation": 189, "longitude": 3.81194}, {"latitude": 47.3553, "elevation": 588, "longitude": 4.77583}, {"latitude": 44.3231, "elevation": 293, "longitude": 4.76222}, {"latitude": 44.8314, "elevation": 46, "longitude": -0.69167}, {"latitude": 47.0586, "elevation": 161, "longitude": 2.35944}, {"latitude": 48.9272, "elevation": 152, "longitude": -0.14944}, {"latitude": 46.6986, "elevation": 157, "longitude": 0.06556}, {"latitude": 43.2167, "elevation": 616, "longitude": 6.37278}, {"latitude": 45.1044, "elevation": 328, "longitude": 1.36944}, {"latitude": 45.29, "elevation": 1115, "longitude": 3.70944}, {"latitude": 43.9906, "elevation": 664, "longitude": 2.60972}, {"latitude": 43.6247, "elevation": 127, "longitude": -0.60917}, {"latitude": 47.3686, "elevation": 904, "longitude": 7.01917}, {"latitude": 48.7158, "elevation": 282, "longitude": 6.58167}, {"latitude": 43.8061, "elevation": 64, "longitude": 4.50278}, {"latitude": 46.0678, "elevation": 899, "longitude": 4.44528}, {"latitude": 42.9183, "elevation": 688, "longitude": 2.865}, {"latitude": 48.4608, "elevation": 96, "longitude": -4.43}, {"latitude": 43.5744, "elevation": 163, "longitude": 1.37611}, {"latitude": 48.7739, "elevation": 168, "longitude": 2.0075}, {"latitude": 47.3375, "elevation": 66, "longitude": -1.65639}, {"latitude": 48.4622, "elevation": 151, "longitude": 4.30944}, {"latitude": 45.8835, "elevation": 246, "longitude": 17.2009}, {"latitude": 45.5027, "elevation": 87, "longitude": 18.5613}, {"latitude": 52.9528, "elevation": 6, "longitude": 4.79061}, {"latitude": 51.8369, "elevation": 2, "longitude": 5.1381}, {"latitude": 50.3942, "elevation": 393, "longitude": 20.0797}, {"latitude": 54.3843, "elevation": 135, "longitude": 18.4563}, {"latitude": 52.4052, "elevation": 92, "longitude": 20.9609}, {"latitude": 50.892, "elevation": 662, "longitude": 16.0395}, {"latitude": 52.4133, "elevation": 94, "longitude": 16.7971}, {"latitude": 50.1517, "elevation": 322, "longitude": 18.7267}, {"latitude": 50.1141, "elevation": 207, "longitude": 22.037}, {"latitude": 53.7903, "elevation": 113, "longitude": 15.8311}, {"latitude": 56.3675, "elevation": 190, "longitude": 12.8517}, {"latitude": 59.6544, "elevation": 50, "longitude": 17.9463}, {"latitude": 57.3034, "elevation": 58, "longitude": 18.4003}, {"latitude": 61.5771, "elevation": 0, "longitude": 16.7144}, {"latitude": 67.7088, "elevation": 0, "longitude": 20.6178}, {"latitude": 56.2955, "elevation": 106, "longitude": 15.6103}, {"latitude": 60.723, "elevation": 0, "longitude": 14.8776}, {"latitude": 65.4309, "elevation": 0, "longitude": 21.865}, {"latitude": 63.6395, "elevation": 64, "longitude": 18.4019}, {"latitude": 63.295, "elevation": 251, "longitude": 14.7591}, {"latitude": 58.2556, "elevation": 136, "longitude": 12.826}, {"latitude": 58.1059, "elevation": 206, "longitude": 15.9363}, {"latitude": 46.0678, "elevation": 937, "longitude": 15.2849}, {"latitude": 46.098, "elevation": 1004, "longitude": 14.2282}, {"latitude": 48.2561, "elevation": 577, "longitude": 17.1531}, {"latitude": 48.7829, "elevation": 1236, "longitude": 20.9873}, {"latitude": 49.2717, "elevation": 1387, "longitude": 19.2493}, {"latitude": 48.2404, "elevation": 613, "longitude": 19.2574}]}');
[tmp.results.elevation]';
%% Export for online visualization
for i_d=1:numel(dc)
    fileID = fopen(['data/exportdDens_' dc(i_d).name '.json'],'w');
    a=struct();
    a.x=datestr(dc(i_d).time,'yyyy-mm-dd HH:MM');
    a.y=dc(i_d).interval*(1/2:double(dc(i_d).levels));
    b = log10(dc(i_d).dens);
    b(dc(i_d).dens==0)=0;
    b(b<0)=0;
    a.data=round(b,2)';
    fprintf(fileID,jsonencode(a));
    fclose(fileID);
end
fileID = fopen('data/exportdRadars.json','w');
fprintf(fileID,jsonencode({dc.name}));
fclose(fileID);


%% Test DBZH -> dens
i_d=14; % seems like a good case study

% eta-> dens
figure; hold on; 
plot([0 max(dc(i_d).dens(:))],[0 max(dc(i_d).dens(:))],'--r');
plot(dc(i_d).eta(:)/11,dc(i_d).dens(:),'.k')

sum(~isnan(dc(i_d).eta(isnan(dc(i_d).dens(:))))) % NO NaN in dens with value in eta
mean(isnan(dc(i_d).eta(dc(i_d).dens(:)==0))) % very few dens=0 -> eta=nan
mean(dc(i_d).eta(dc(i_d).dens(:)==0)==0) % very few dens=0 -> eta=0

id = dc(i_d).sd_vvp(:)>2;
figure; 
plot(dc(i_d).DBZH(id),dc(i_d).sd_vvp(id),'.k')
legend(num2str(sum(id)))

scatter(dc(i_d).eta(:)/11,dc(i_d).dens(:),[],dc(i_d).sd_vvp(:)<2,'.')


% DBZH -> eta
figure; hold on;
a=361*10.^(dc(i_d).DBZH(:)/10);
b=dc(i_d).eta(:);
plot([0 max(dc(i_d).eta(:))],[0 max(dc(i_d).eta(:))],'--r');
plot(a,b,'.k');
axis([0 max(dc(i_d).eta(:)) 0 max(dc(i_d).eta(:))])

% Rain between DBZH and eta
figure; plot(dc(i_d).eta(:),dc(i_d).DBZH(:),'.')

figure; hold on;
histogram(dc(i_d).DBZH(isnan(dc(i_d).eta(:))))
histogram(dc(i_d).DBZH(~isnan(dc(i_d).eta(:))))
