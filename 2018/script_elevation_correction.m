% script


load('data\dc_corr'); addpath('functions'); load coastlines.mat

%% Analyse of the vertical profile
fig=figure; hold on;
isdecreasing=false(numel(dc),2);
for i_d=1:numel(dc)
    id = dc(i_d).scatter_lim:25;
    tmp = dc(i_d).dens4(:,id);%./repmat(nansum(dc(i_d).dens4(:,id),2),1,numel(id));
    
    subplot(1,2,1); hold on;
    it1=dc(i_d).time>datetime('01-March-2018') & dc(i_d).time<datetime('01-May-2018');
    plot(nanmean(tmp(it1,:)),dc(1).alt(id)-dc(i_d).heightDEM);
    axis([0 30 0 5000]); title('Spring'); xlabel('Profile % of total dens'); ylabel('altitude [m]')
    
    subplot(1,2,2); hold on;
    it2=dc(i_d).time>datetime('01-Sep-2018') & dc(i_d).time<datetime('01-Nov-2018');
    plot(nanmean(tmp(it2,:)),dc(1).alt(id)-dc(i_d).heightDEM);
    axis([0 30 0 5000]); title('Autumn');xlabel('Profile % of total dens'); ylabel('altitude [m]')
    
    %       plot(nanmean(tmp),dc(1).alt(id)-dc(i_d).heightDEM);
    %       axis([0 50 0 5000]); title('Year round'); xlabel('Profile % of total dens'); ylabel('altitude [m]')
    
end
legend({dc.name})

datacursormode on;
dcm = datacursormode(fig);
c_info = getCursorInfo(dcm);
% disp(get(c_info.Target, 'DisplayName'));
set(dcm,'UpdateFcn',@customdatatip)

figure;
worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow('worldrivers.shp','Color', 'blue')
scatterm([dc.lat],[dc.lon],100,[dc.heightDEM],'filled','MarkerEdgeColor','k');

figure; hold on;
plot(coastlon,coastlat,'k')
for i_d=1:numel(dc)
    id = dc(i_d).scatter_lim:25;
    tmp = dc(i_d).dens4(:,id);%./repmat(nansum(dc(i_d).dens4(:,id),2),1,numel(id));

    it1=dc(i_d).time>datetime('01-March-2018') & dc(i_d).time<datetime('01-May-2018');
    scatter3(dc(i_d).lon*ones(numel(id),1), dc(i_d).lat*ones(numel(id),1), dc(1).alt(id)-dc(i_d).heightDEM, [], nanmean(tmp(it1,:)), 'filled');
end
xlim([min([dc.lon]) max([dc.lon])]); ylim([min([dc.lat]) max([dc.lat])]); zlim([0 2000])


%% Extrapolate elevation with fitted function: NOT USED
ft = fittype('a1*exp(-b*x)-a2*exp(-b*x)');

x=0:10:5000;
a=50;ar=0.3;b1=1000;b2=100;
y=a*exp(-x/b1)-a*ar*exp(-x/b2);
plot(y,x);

for i_d=1:numel(dc)
    id = dc(i_d).scatter_lim:25;
    it1=dc(i_d).time>datetime('01-March-2018') & dc(i_d).time<datetime('01-May-2018') & ~all(isnan(dc(i_d).dens4(:,id)),2)';
    tmp = dc(i_d).dens4(it1,id);
    alt = dc(i_d).alt(id)-dc(i_d).heightDEM;
    
    for i_alt=1:numel(alt)
        subplot(5,5,i_alt)
        ksdensity(tmp(:,i_alt),linspace(0,50,100))
    end
    
    
    f = @(A,alt,ar,b1,b2) bsxfun(@times,A(:)',exp(-alt(:)./b1))-bsxfun(@times,A(:)',ar.*exp(-alt(:)./b2));
    
    A = nanmax(tmp,[],2);
    f_min_1= @(x)sqrt(nansum(( nanmean(tmp./repmat(A,1,numel(id)))' - f(x(4:end), alt ,x(1),x(2),x(3)) ).^2 ));
    x0tmp=fminsearch(f_min_1 ,[0 1000 100 1]);
    
    figure; hold on;
    plot(nanmean(tmp./repmat(A,1,numel(id)))',alt);
    plot(f(x0tmp(4:end), 0:10:5000 ,x0tmp(1),x0tmp(2),x0tmp(3)), 0:10:5000);
    
    
    f_min_2 = @(x) sqrt(nansum(nansum((tmp' - f(x, alt , x0tmp(1), x0tmp(2), x0tmp(3)) ).^2 )));
    x=fminsearch(f_min_2,x0tmp(end).*A(:)',optimset('MaxFunEvals',1000));
    

end


%% MPS filling
% simulate the missing data below scatter_lim of each radar. Simulation is
% per week using a training image of  +/- 2 weeks.


serverAddress='wally-front1.unil.ch'; %wally-front1.unil.ch tesla-k20c.gaia.unil.ch knl.gaia.unil.ch
ti_size=500; % number of nights.radar to consider (sampled randomly among all available
ki_size=[49 7]; % windows of search, (1) grid size is 5min, (2) 200m
n_real=30;
max_alt = 5; % 5-> 1000m at 200 resolution
max_dist = 600; %km only select nearby radar 
MPS_full=cell(size(dc));

Ddist=squareform(pdist([[dc.lat]' [dc.lon]'],@lldistkm));

scoret = datenum(repmat(dc(1).time',1,numel(dc))-mean(cat(3,[dc.sunrise],[dc.sunset]),3)) ./ datenum([dc.sunrise]-[dc.sunset])*2;
tti = cell(1,numel(dc));
ddi = cell(1,numel(dc));
for i_d=1:numel(dc)
    ddi{i_d} = dc(i_d).dens4;
    ddi{i_d}(:,1:dc(i_d).scatter_lim-1)=nan;
    
    id = dc(i_d).alt > dc(i_d).heightDEM;
    tti{i_d} = ddi{i_d}(:,id);
    
    MPS_full{i_d}=nan(numel(dc(i_d).time),size(dc(i_d).dens4,2),n_real);
end

days7 = datetime('31-December-2017 12:00'):7:datetime('01-January-2019 12:00');


for i_t=1:numel(days7)-1

    % Find the corresponding time
    i_t_di = dc(1).time>days7(i_t) & dc(1).time<days7(i_t+1);
    % Select training set
    ti = cellfun(@(x) log10(x(i_t_di,1:max_alt)+1),tti,'UniformOutput',false);
    
    for i_d= 1:numel(dc)
        di = log10(ddi{i_d}(i_t_di,:)+1);
        
        % Define path
        % sp = reshape((size(di,1)*size(di,2)):-1:1,size(di,1),[]);
        sp = repmat(size(di,2):-1:1,size(di,1),1);
        sp(all(isnan(di),2),:)=-Inf;
        sp(:,dc(i_d).scatter_lim:end)=-Inf;
        
        % Subsample nearby radar
        tii = ti(Ddist(i_d,:)>0 & Ddist(i_d,:)<max_dist);
        
        for i_r=1:n_real
            if any(~isinf(sp(:)))
  
                sp2 = sp+rand(size(sp))*0.9;
                
                % [result,t]=g2s('-sa','knl.gaia.unil.ch','-a','qs','-ti',ti{unique(unidrnd(numel(ti),ti_size,1))},'-di',di,'-dt',zeros(1,1),'-k',1.5,'-n',15,'-sp',sp,'-ki',ones(ki_size(1),ki_size(2)),'-j',1,64,1);
                % MPS{i_d}(i_t_di,:,i_real) = result(:,1:dc(i_d).scatter_lim);

    %             [result,t]=g2s('-sa','tesla-k20c.gaia.unil.ch','-a','qs','-ti',ti,'-di',di,'-dt',zeros(1,1),'-k',1.5,'-n',15,'-sp',sp,'-ki',ones(ki_size(1),ki_size(2)),'-j',1,20,1);
    %              MPS{i_d}(i_t_di,:,:) = result(:,1:dc(i_d).scatter_lim,:);

                % id(i_t,i_d)=g2s('-sa','wally-front1.unil.ch','-a','qs','-ti',ti{:},'-di',di,'-dt',zeros(1,1),'-k',1.5,'-n',15,'-sp',sp,'-ki',ones(ki_size(1),ki_size(2)),'-j',1,16,1,'-submitOnly','-account','gmariet1_gaia','-after',0);
                % result=g2s('-sa','wally-front1.unil.ch','-waitAndDownload',id(i_t,i_d));
                % MPS{i_d}(i_t_di,:,:) = result(:,1:dc(i_d).scatter_lim,:);

                MPS_full{i_d}(i_t_di,:,i_r) = g2s('-sa','knl.gaia.unil.ch','-a','qs','-ti',tii{:},'-di',di,'-dt',zeros(1,1),'-k',1.5,'-n',15,'-sp',sp2,'-ki',ones(ki_size(1),ki_size(2)),'-j',8,8,1);
                % MPS{i_d}(i_t_di,:,i_r) = 10.^(result(:,1:dc(i_d).scatter_lim-1))-1;
            else
                MPS_full{i_d}(i_t_di,:,i_r)=nan;
            end
        end
    end
end

for i_d=1:numel(dc)
    MPS{i_d} = 10.^( MPS_full{i_d}(:,1:dc(i_d).scatter_lim,:)) -1;
    if (size(MPS{i_d},2)==1)
        MPS{i_d}(:,1,:)=[];
    end
end


% figure test to check everything is ok. 
figure;
subplot(3,1,1); imagesc(tii{3}'); set(gca,'ydir','normal')
subplot(3,1,2); imagesc(di'); set(gca,'ydir','normal')
% subplot(3,1,1); imagesc(sp'); set(gca,'ydir','normal')
subplot(3,1,3); imagesc(result'); set(gca,'ydir','normal')

% Send
save('data/BelowRadarMPS4','MPS','-v7.3');
% MPS, MPS and MPS can be then deleted. 


load('./data/BelowRadarMPS.mat'); 

fig=figure; hold on;
it1=dc(1).time>datetime('01-March-2018') & dc(1).time<datetime('01-May-2018');
it2=dc(1).time>datetime('01-Sep-2018') & dc(1).time<datetime('01-Nov-2018');
for i_d=1:numel(dc)
    %clf;
    if size(MPS{i_d},2)~= 0
        subplot(1,2,1); hold on;
        plot(nanmean(dc(i_d).dens4(it1,:)),dc(1).alt);
        plot(nanmean(mean(MPS{i_d}(it1,:,:),3)),dc(1).alt(1:dc(i_d).scatter_lim),'--');
        plot(nanmean(quantile(MPS{i_d}(it1,:,:),0.2,3)),dc(1).alt(1:dc(i_d).scatter_lim),'--');
        plot(nanmean(quantile(MPS{i_d}(it1,:,:),0.8,3)),dc(1).alt(1:dc(i_d).scatter_lim),'--');
        axis([0 30 0 5000]); title('Spring'); xlabel('Profile  dens'); ylabel('altitude [m]')

        subplot(1,2,2); hold on;
        plot(nanmean(dc(i_d).dens4(it2,:)),dc(1).alt);
        plot(nanmean(mean(MPS{i_d}(it2,:,:),3)),dc(1).alt(1:dc(i_d).scatter_lim),'--');
        plot(nanmean(quantile(MPS{i_d}(it2,:,:),0.2,3)),dc(1).alt(1:dc(i_d).scatter_lim),'--');
        plot(nanmean(quantile(MPS{i_d}(it2,:,:),0.8,3)),dc(1).alt(1:dc(i_d).scatter_lim),'--');
        axis([0 30 0 5000]); title('Autumn');xlabel('Profile  dens'); ylabel('altitude [m]')
    end
    % keyboard;
end


fig=figure; hold on;
it(1,:)=dc(1).time>datetime('01-March-2018') & dc(1).time<datetime('01-May-2018');
it(2,:)=dc(1).time>datetime('01-Sep-2018') & dc(1).time<datetime('01-Nov-2018');
for i_d=1:numel(dc)
    %clf;
    if size(MPS{i_d},2)~= 0
        for s=1:2
            subplot(1,2,s); hold on;
            plot(nanmean(dc(i_d).dens4(it(s,:),:)),dc(1).alt-dc(i_d).heightDEM);
            y = dc(1).alt(1:dc(i_d).scatter_lim)-dc(i_d).heightDEM;
            x1=nanmean(quantile(MPS{i_d}(it(s,:),:,:),0.1,3));
            x2=nanmean(quantile(MPS{i_d}(it(s,:),:,:),0.9,3));
            fill( [x1 fliplr(x2)], [y fliplr(y)],'k')
            plot(nanmean(mean(MPS{i_d}(it(s,:),:,:),3)),y,'--');
            axis([0 30 -500 2000]); title('Spring'); xlabel('Profile  dens'); ylabel('altitude [m]')
        end
    end
end


%% Get elevatation around

% Find elevation with API OPEN_ELEVATION. NOT ALWAYS CORRECT. EXECEL
% spreadsheet has corrected value with https://www.freemaptools.com/elevation-finder.htm
% strloc = join(strcat( strcat(cellstr(num2str([dc.lat]')),',') , cellstr(num2str([dc.lon]')) )','|');
% strloc = replace(strloc{1},' ','');
% tmp = webread(['https://api.open-elevation.com/api/v1/lookup?locations=' strloc ]); % num2str(d(i_d).lat) ',' num2str(d(i_d).lon)
% tmp = jsondecode('{"results": [{"latitude": 51.1917, "elevation": -1, "longitude": 3.0642}, {"latitude": 49.9143, "elevation": 546, "longitude": 5.5056}, {"latitude": 49.6583, "elevation": 861, "longitude": 13.8178}, {"latitude": 49.5011, "elevation": 744, "longitude": 16.7885}, {"latitude": 53.564, "elevation": 1, "longitude": 6.7483}, {"latitude": 54.0044, "elevation": 93, "longitude": 10.0468}, {"latitude": 51.1246, "elevation": 225, "longitude": 13.7686}, {"latitude": 49.5407, "elevation": 774, "longitude": 12.4028}, {"latitude": 53.3394, "elevation": 4, "longitude": 7.025}, {"latitude": 51.4055, "elevation": 153, "longitude": 6.9669}, {"latitude": 51.3112, "elevation": 560, "longitude": 8.802}, {"latitude": 47.8736, "elevation": 1487, "longitude": 8.0036}, {"latitude": 52.4601, "elevation": 58, "longitude": 9.6945}, {"latitude": 48.0431, "elevation": 668, "longitude": 10.2204}, {"latitude": 50.5001, "elevation": 844, "longitude": 11.135}, {"latitude": 50.1097, "elevation": 553, "longitude": 6.5485}, {"latitude": 49.9859, "elevation": 204, "longitude": 8.714}, {"latitude": 52.6486, "elevation": 145, "longitude": 13.8578}, {"latitude": 54.1757, "elevation": 2, "longitude": 12.0581}, {"latitude": 48.1747, "elevation": 633, "longitude": 12.1018}, {"latitude": 48.5853, "elevation": 733, "longitude": 9.7828}, {"latitude": 52.1601, "elevation": 153, "longitude": 11.1761}, {"latitude": 55.1127, "elevation": 162, "longitude": 14.8875}, {"latitude": 55.1731, "elevation": 2, "longitude": 8.552}, {"latitude": 57.4893, "elevation": 85, "longitude": 10.1365}, {"latitude": 55.3262, "elevation": 31, "longitude": 12.4493}, {"latitude": 56.024, "elevation": 120, "longitude": 10.0246}, {"latitude": 59.3977, "elevation": 33, "longitude": 24.6021}, {"latitude": 58.4823, "elevation": 131, "longitude": 25.5187}, {"latitude": 36.8325, "elevation": 479, "longitude": -2.08222}, {"latitude": 39.4289, "elevation": 605, "longitude": -6.28528}, {"latitude": 41.4081, "elevation": 613, "longitude": 1.88472}, {"latitude": 43.1689, "elevation": 592, "longitude": -8.52694}, {"latitude": 41.9956, "elevation": 873, "longitude": -4.60278}, {"latitude": 28.0186, "elevation": 1760, "longitude": -15.6144}, {"latitude": 40.1758, "elevation": 702, "longitude": -3.71361}, {"latitude": 36.6133, "elevation": 1131, "longitude": -4.65917}, {"latitude": 38.2644, "elevation": 1187, "longitude": -1.18972}, {"latitude": 39.3797, "elevation": 109, "longitude": 2.785}, {"latitude": 43.4625, "elevation": 867, "longitude": -6.30194}, {"latitude": 37.6875, "elevation": 501, "longitude": -6.33444}, {"latitude": 43.4033, "elevation": 554, "longitude": -2.84194}, {"latitude": 39.1761, "elevation": 166, "longitude": -0.25211}, {"latitude": 41.7339, "elevation": 811, "longitude": -0.54583}, {"latitude": 60.9039, "elevation": 181, "longitude": 27.1081}, {"latitude": 61.7673, "elevation": 157, "longitude": 23.0764}, {"latitude": 61.907, "elevation": 186, "longitude": 29.7977}, {"latitude": 60.1285, "elevation": 113, "longitude": 21.6434}, {"latitude": 62.8626, "elevation": 107, "longitude": 27.3815}, {"latitude": 67.1391, "elevation": 0, "longitude": 26.8969}, {"latitude": 64.7749, "elevation": 93, "longitude": 26.3189}, {"latitude": 60.2706, "elevation": 160, "longitude": 24.869}, {"latitude": 63.1048, "elevation": 0, "longitude": 23.8209}, {"latitude": 50.1358, "elevation": 69, "longitude": 1.83472}, {"latitude": 42.1297, "elevation": 40, "longitude": 9.49639}, {"latitude": 50.1283, "elevation": 189, "longitude": 3.81194}, {"latitude": 47.3553, "elevation": 588, "longitude": 4.77583}, {"latitude": 44.3231, "elevation": 293, "longitude": 4.76222}, {"latitude": 44.8314, "elevation": 46, "longitude": -0.69167}, {"latitude": 47.0586, "elevation": 161, "longitude": 2.35944}, {"latitude": 48.9272, "elevation": 152, "longitude": -0.14944}, {"latitude": 46.6986, "elevation": 157, "longitude": 0.06556}, {"latitude": 43.2167, "elevation": 616, "longitude": 6.37278}, {"latitude": 45.1044, "elevation": 328, "longitude": 1.36944}, {"latitude": 45.29, "elevation": 1115, "longitude": 3.70944}, {"latitude": 43.9906, "elevation": 664, "longitude": 2.60972}, {"latitude": 43.6247, "elevation": 127, "longitude": -0.60917}, {"latitude": 47.3686, "elevation": 904, "longitude": 7.01917}, {"latitude": 48.7158, "elevation": 282, "longitude": 6.58167}, {"latitude": 43.8061, "elevation": 64, "longitude": 4.50278}, {"latitude": 46.0678, "elevation": 899, "longitude": 4.44528}, {"latitude": 42.9183, "elevation": 688, "longitude": 2.865}, {"latitude": 48.4608, "elevation": 96, "longitude": -4.43}, {"latitude": 43.5744, "elevation": 163, "longitude": 1.37611}, {"latitude": 48.7739, "elevation": 168, "longitude": 2.0075}, {"latitude": 47.3375, "elevation": 66, "longitude": -1.65639}, {"latitude": 48.4622, "elevation": 151, "longitude": 4.30944}, {"latitude": 45.8835, "elevation": 246, "longitude": 17.2009}, {"latitude": 45.5027, "elevation": 87, "longitude": 18.5613}, {"latitude": 52.9528, "elevation": 6, "longitude": 4.79061}, {"latitude": 51.8369, "elevation": 2, "longitude": 5.1381}, {"latitude": 50.3942, "elevation": 393, "longitude": 20.0797}, {"latitude": 54.3843, "elevation": 135, "longitude": 18.4563}, {"latitude": 52.4052, "elevation": 92, "longitude": 20.9609}, {"latitude": 50.892, "elevation": 662, "longitude": 16.0395}, {"latitude": 52.4133, "elevation": 94, "longitude": 16.7971}, {"latitude": 50.1517, "elevation": 322, "longitude": 18.7267}, {"latitude": 50.1141, "elevation": 207, "longitude": 22.037}, {"latitude": 53.7903, "elevation": 113, "longitude": 15.8311}, {"latitude": 56.3675, "elevation": 190, "longitude": 12.8517}, {"latitude": 59.6544, "elevation": 50, "longitude": 17.9463}, {"latitude": 57.3034, "elevation": 58, "longitude": 18.4003}, {"latitude": 61.5771, "elevation": 0, "longitude": 16.7144}, {"latitude": 67.7088, "elevation": 0, "longitude": 20.6178}, {"latitude": 56.2955, "elevation": 106, "longitude": 15.6103}, {"latitude": 60.723, "elevation": 0, "longitude": 14.8776}, {"latitude": 65.4309, "elevation": 0, "longitude": 21.865}, {"latitude": 63.6395, "elevation": 64, "longitude": 18.4019}, {"latitude": 63.295, "elevation": 251, "longitude": 14.7591}, {"latitude": 58.2556, "elevation": 136, "longitude": 12.826}, {"latitude": 58.1059, "elevation": 206, "longitude": 15.9363}, {"latitude": 46.0678, "elevation": 937, "longitude": 15.2849}, {"latitude": 46.098, "elevation": 1004, "longitude": 14.2282}, {"latitude": 48.2561, "elevation": 577, "longitude": 17.1531}, {"latitude": 48.7829, "elevation": 1236, "longitude": 20.9873}, {"latitude": 49.2717, "elevation": 1387, "longitude": 19.2493}, {"latitude": 48.2404, "elevation": 613, "longitude": 19.2574}]}');
% [tmp.results.elevation]';

files={'gt30w020n90'};%,'gt30e020n90'}; % ,'gt30e020n40','gt30w020n40'
[DEM,GeoCellRef] = geotiffread(['data/' files{1}]);
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
V = nan(numel(alpha),numel(dc(1).alt),numel(dc));

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
        % elevation in the direction of alpha
        DEM_id_r = imrotate(DEM_id,alpha(i_alpha));
        % We take the max altitude (DEM_id_r) on the path perpendicular of the transect.
        % biorad compute the density over the available volume (e.g. remove
        % mountain). But the volume available for flight is still
        % different.
        tmp = max(double(DEM_id_r)); % max(max(double(DEM_id_r)),dc(i_d).height);
        % Compute the the altitude of available flight (distance from ground to top of bin height).
        tmp2 = max(min(bsxfun(@minus,dc(1).alt'+dc(i_d).interval/2,tmp),dc(i_d).interval),0)./1000;
        
        % The area is not homogenously sampled along the transect as it is a circle. We normalized
        % it with J
        coeff_circle = sum(J)/sum(J(:));
        V(i_alpha,:,i_d) = sum( tmp2 .* repmat(coeff_circle,numel(dc(1).alt),1) ,2);
        
        subplot(5,4,i_alpha); hold on;
        % imagesc(DEM_id_r);
        surf(1:size(tmp2,2),dc(i_d).alt,tmp2)
        plot(tmp,'linewidth',3); %plot(0,[0 dc(i_d).alt+100],'xr')
        title(['a=' num2str(alpha(i_alpha)) '   V=' num2str(sum(V(i_alpha,:,i_d)))])
    end
    dc(i_d).VolBelow = V(:,:,i_d);
end

% Figure

dc(i_d).u = dc(i_d).ff .* sind(dc(i_d).dd); % m/s | 0° is north and 90° is west. -> u is east (+) - west (-)
dc(i_d).v = dc(i_d).ff .* cosd(dc(i_d).dd);

figure;
subplot(1,4,[1 3]);hold on; title(dc(i_d).name)
%scatter3(GeoCellRef_LAT(:), GeoCellRef_LON(:), DEM(:),'.k')
surf(GeoCellRef_LON(id_any_1, id_any_2), GeoCellRef_LAT(id_any_1, id_any_2), DEM_id,'EdgeColor','none')
plot3(dc(i_d).lon, dc(i_d).lat, dc(i_d).height,'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10);
[X,Y,Z] = cylinder(dc(i_d).maxrange/111);
Z(2,:) = dc(i_d).height;
mesh(X+dc(i_d).lon,Y+dc(i_d).lat,Z,'FaceAlpha',0)
xlabel('lat'); ylabel('lon'); zlabel('elev');
if all(DEM_id(:)==0)
    zlim0=0;
else
    zlim0=min(DEM_id(~DEM_id==0));
end
zlim0=min(zlim0,dc(i_d).height);
zlim([zlim0 dc(i_d).height+700]);axis_t=axis; view(0,0)
id = dc(i_d).scatter_lim:25;
it1=dc(i_d).time>datetime('01-March-2018') & dc(i_d).time<datetime('01-May-2018');
it2=dc(i_d).time>datetime('01-Sep-2018') & dc(i_d).time<datetime('01-Nov-2018');
mu=nanmean(dc(i_d).u(it1,id) .* dc(i_d).dens4(it1,id) ./ repmat(nanmean(dc(i_d).dens4(it1,id)),sum(it1),1) );
mv=nanmean(dc(i_d).v(it1,id) .* dc(i_d).dens4(it1,id) ./ repmat(nanmean(dc(i_d).dens4(it1,id)),sum(it1),1) );
a = [zeros(numel(id),1) mu'/20 nan(numel(id),1)]';
b = [zeros(numel(id),1) mv'/20 nan(numel(id),1)]';
c = [dc(1).alt(id)' dc(1).alt(id)' nan(numel(id),1)]';
plot3(dc(i_d).lon+a(:), dc(i_d).lat+b(:), c(:),'g','linewidth',3)
mu=nanmean(dc(i_d).u(it2,id) .* dc(i_d).dens4(it2,id) ./ repmat(nanmean(dc(i_d).dens4(it2,id)),sum(it2),1) );
mv=nanmean(dc(i_d).v(it2,id) .* dc(i_d).dens4(it2,id) ./ repmat(nanmean(dc(i_d).dens4(it2,id)),sum(it2),1) );
a = [zeros(numel(id),1) mu'/1000*60*60/72 nan(numel(id),1)]';
b = [zeros(numel(id),1) mv'/1000*60*60/111 nan(numel(id),1)]';
plot3(dc(i_d).lon+a(:), dc(i_d).lat+b(:), c(:),'r','linewidth',3)
axis(axis_t); set(gca,'DataAspectRatio',[1 1 300])
subplot(1,4,4);  hold on;
plot(nanmean(dc(i_d).dens4(it1,id)),dc(1).alt(id),'g','linewidth',3)
plot(nanmean(dc(i_d).dens4(it2,id)),dc(1).alt(id),'r','linewidth',3);
plot(nanmean(dc(i_d).dens4(:,id)),dc(1).alt(id));
plot(0, dc(i_d).height,'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10);
ylim([zlim0 dc(i_d).height+700]);xlim([0 25])


% Figure volume as a function of direction. 
figure;
for i_d=1:numel(dc)
    vol_old =sum(diff(max( [0; dc(1).alt'+dc(i_d).interval/2]-dc(i_d).height,0)));
    subplot(6,7,i_d);
    x1=deg2rad(alpha); y1=sum(dc(i_d).VolBelow,2)*1000-vol_old; x1(y1<=0)=NaN; y1(y1<=0)=NaN;
    polarplot(x1,abs(y1))
    hold on;
    x1=deg2rad(alpha); y1=sum(dc(i_d).VolBelow,2)*1000-vol_old; x1(y1>=0)=NaN; y1(y1>=0)=NaN;
    polarplot(x1,abs(y1));
    ax=gca;
    ax.ThetaDir = 'clockwise';
    ax.ThetaLim=[-90 90];
    ax.ThetaZeroLocation='top';
    % rlim([0 1.2])
    %thetaticks('');
    title(dc(i_d).name)
end


%Save
for i_d=1:numel(dc)
    VolBelow(i_d,:,:) = dc(i_d).VolBelow;
end
save('data/VolBelow','VolBelow','alpha');
save('data/dc_corr','dc','start_date','end_date','quantity','-v7.3')


%% Figures
load('data/dc_corr')
load('data/VolBelow')
load('./data/BelowRadarMPS.mat'); 

% Figure: altitudinal coverage of the data
V=nan(25,37);
for i_d=1:numel(dc)
    V(:,i_d)=mean(dc(i_d).VolBelow)/.2;
end
V=1-V;
figure; hold on;
surf(1:37,dc(1).alt,V)
tmp = linspace(1,0,64)';
colormap( tmp.*ones(64,3) + (1-tmp).*flipud(summer))%demcmap()
plot(1:37,[dc.height],'ok','MarkerFaceColor','k')
stairs(0.5:37.5,[dc.heightDEM dc(end).heightDEM],'-k','linewidth',3)
stairs(0.5:37.5,dc(1).alt([dc.scatter_lim dc(end).scatter_lim])-100,'--r','linewidth',2)
xticks(1:37)
xticklabels({dc.name})
yticks(dc(1).alt)
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';
xtickangle(90); axis tight
ylim([0 1500])
legend({'Radar','Ground','First bin with acceptable data'})
c=colorbar; c.Label.String='Proportion of alt. bin under ground';




% Estimate total influence of elevation correction on MTR
% Initiate empty density matrix
data_denss_sim = nan(numel(dc(1).time),  numel(dc), size(MPS{1},3));
data_MTRs_sim = data_denss_sim;
data_denss = nan(numel(dc(1).time),  numel(dc));
data_MTRs = data_denss;

% Loop for eac radar
for i_d=1:numel(dc)
    % Vertical agregation accounts for the DEM based on VolBelow.
    % First, find the direction of flight of the lowest non-nan bin
    [~, id]  = max(~isnan(dc(i_d).dd), [], 2);
    dir_lowest = dc(i_d).dd(sub2ind(size(dc(i_d).dd), (1:numel(dc(1).time))',id));
    % Transform this direction angle (0-360°) into the index of alpha (-90°-90°)
    dir_lowest_round = round(dir_lowest,-1);
    id_alpha=nan(size(dir_lowest_round));
    id_alpha(dir_lowest_round>=0 & dir_lowest_round<90) = dir_lowest_round(dir_lowest_round>=0 & dir_lowest_round<90)+90;
    id_alpha(dir_lowest_round>=90 & dir_lowest_round<270) = dir_lowest_round(dir_lowest_round>=90 & dir_lowest_round<270)-90;
    id_alpha(dir_lowest_round>=270 & dir_lowest_round<360) = dir_lowest_round(dir_lowest_round>=270 & dir_lowest_round<360)-270;
    id_alpha = id_alpha/10+1;
    height_vol = nan(numel(dc(1).time),numel(dc(1).alt));
    height_vol(~isnan(id_alpha),:) = dc(i_d).VolBelow(id_alpha(~isnan(id_alpha)),:);
    
    % Integrate volume density into surface density based on the direction and corresponding volume. bird/km^2
    % data.denss(:,i_d) = sum(data_dens(:,:,i_d) .* height_vol,2);

   
    % Interpolation speed down to the ground
    v_g = sqrt( dc(i_d).vb.^2 + dc(i_d).ub.^2 );
    id_min = find(diff(all(isnan(v_g)))==-1);
    v_g(:,1:id_min)= repmat(v_g(:,id_min+1),1,id_min);
    
    % Computing with result of MPS
    tmp = dc(i_d).dens4;
    data_denss_MPS = nan(size(tmp,1),size(MPS{i_d},3),size(tmp,2));
    for i_real=1:size(MPS{i_d},3)
        tmp(:,1:dc(i_d).scatter_lim-1) = MPS{i_d}(:,1:dc(i_d).scatter_lim-1,i_real);
        data_denss_MPS(:,i_real,:) = tmp .* height_vol;
        data_MTRs_sim(:,i_d,i_real) = nansum(tmp .* height_vol .* v_g,2);
    end
    data_denss_sim(:,i_d,:) = sum(data_denss_MPS,3);
     
    % Standar method
    vidh = min(max(dc(1).alt+100-dc(i_d).heightDEM,0),200)/1000;
    data_denss(:,i_d) = sum(dc(i_d).dens4 .* repmat(vidh,numel(dc(i_d).time),1),2);
    data_MTRs(:,i_d) = nansum(dc(i_d).dens4 .* repmat(vidh,numel(dc(i_d).time),1) .* v_g,2);
end

diff_dens = data_denss_sim - repmat(data_denss,1,1,size(MPS{i_d},3));
diff_MTR = data_MTRs_sim - repmat(data_MTRs,1,1,size(MPS{i_d},3));

% Figure
tmp2s = nan(3,numel(dc));
tmp2a = nan(3,numel(dc));
figure;
for i_d=1:numel(dc)
    % subplot(6,7,i_d)
    % histogram(nansum(diff_dens(:,i_d,:))./nansum(data_denss(:,i_d)))
    tmps=nansum(diff_MTR(1:end/2,i_d,:))./nansum(data_MTRs(1:end/2,i_d));
    tmpa=nansum(diff_MTR(end/2:end,i_d,:))./nansum(data_MTRs(end/2:end,i_d));
    tmp2s(:,i_d) = [mean(tmps) min(tmps) max(tmps)];
    tmp2a(:,i_d) = [mean(tmpa) min(tmpa) max(tmpa)];
    %histogram(tmp*100)
    %xlabel(dc(i_d).name); xtickformat('percentage')
end

figure; hold on;
plot([0 0],[1 37],'--k')
tmp3=tmp2s*100;
errorbar(tmp3(1,:),1:numel(dc),[],[],tmp3(2,:)-tmp3(1,:),tmp3(3,:)-tmp3(1,:),'.k');
scatter(tmp3(1,:),1:numel(dc),[],'o','filled')
tmp3=tmp2a*100;
errorbar(tmp3(1,:),1:numel(dc),[],[],tmp3(2,:)-tmp3(1,:),tmp3(3,:)-tmp3(1,:),'.k');
scatter(tmp3(1,:),1:numel(dc),[],'o','filled')
yticks(1:37); yticklabels({dc.name}); xtickformat('percentage'); ylim([0 38])
ylabel('Radars'); xlabel('Relative change of mean MTR between new and old method')
grid on; box on










