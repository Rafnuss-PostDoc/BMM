%% Birds at risk
% Determine how many birds are at risk of collision with windturbine where and 
% when.
%% Setup
% Load data

load('data/windfarms_processed.mat')
load('data/energy_processed.mat');
load('../2018/data/Density_estimationMap','g');
load('../2018/data/Density_estimationMap','gd');
load('../2018/data/Flight_estimationMap','guv');
load('./data/ratio');
addpath('functions')

%% 
% Compute area
dy = lldistkm([g.lat(1) g.lon(1)],[g.lat(2) g.lon(1)]);
dx = lldistkm([g.lat2D(:,1) g.lon2D(:,1)],[g.lat2D(:,1) g.lon2D(:,2)]);
area = repmat(dx*dy,1,g.nlon);

% Create country map
bd = load('borderdata.mat');
inCountry=false(numel(g.lat),numel(g.lon),3);
ind = strcmp(bd.places,'France');
inCountry(:,:,1) = reshape(inpolygon(g.lon2D(:), g.lat2D(:), bd.lon{ind}', bd.lat{ind}'),size(g.lon2D));
ind = strcmp(bd.places,'Germany');
inCountry(:,:,2) = reshape(inpolygon(g.lon2D(:), g.lat2D(:), bd.lon{ind}', bd.lat{ind}'),size(g.lon2D));
ind =  strcmp(bd.places,{'Belgium'}) |  strcmp(bd.places,{'Netherlands'}) |  strcmp(bd.places,{'Luxembourg'});
inCountry(:,:,3) = reshape(inpolygon(g.lon2D(:), g.lat2D(:), [bd.lon{ind}]', [bd.lat{ind}]'),size(g.lon2D));

% load coastlines
load coastlines
dcoast = pdist2([g.lat2D(:) g.lon2D(:)],[coastlat coastlon], @lldistkm);
dcoast = reshape(min(dcoast,[],2),size(g.lat2D));

% Define season
season = nan(size(time));
season(time<datetime('1-March-2018')|time>=datetime('1-Dec-2018'))=1;
season(time>=datetime('1-March-2018')&time<datetime('1-May-2018'))=2;
season(time>=datetime('1-May-2018')&time<datetime('1-Sep-2018'))=3;
season(time>=datetime('1-Sep-2018')&time<datetime('1-Dec-2018'))=4;


%% Construct turbine sweep area 
[G,ID] = findgroups(tim.grid_id);
turbineSweptArea=nan(1,g.nlm);
turbineSweptArea(ID) = splitapply(@nansum,tim.turbine_sweptArea.*tim.Number_of_turbines,G);
turbineSweptAreaMap = nan(size(g.latlonmask));
turbineSweptAreaMap(g.latlonmask) = turbineSweptArea;


%% Compute birdDir at 1hr resolution
[~,Locb]=ismember(dateshift(g.time,'end','hour'),time);
birdDir_tmp = nan(g.nlat,g.nlon,g.nt);
birdDir_tmp(repmat(g.latlonmask,1,1,g.nt)) = angle(guv.u_est + 1i*guv.v_est);
birdDir = nan(g.nlat,g.nlon,numel(time));
for i=1:numel(time)
    birdDir(:,:,i) = meanangle(birdDir_tmp(:,:,Locb==i),3);
end

%% Windspeed

file='./data/ECMWF/2018_srf.nc'; %ncdisp(file);
w_time = datetime('2018-01-01 00:00'):1/24:datetime('2018-12-31 23:00'); % ncread(file,'time')
w_lat=flip(double(ncread(file,'latitude')));
w_lon=double(ncread(file,'longitude'));
w_u100 = permute(flip(ncread(file,'u100'),2) , [2 1 3]); 
w_v100 = permute(flip(ncread(file,'v100'),2) , [2 1 3]);
w_u10 = permute(flip(ncread(file,'u10'),2) , [2 1 3]); 
w_v10 = permute(flip(ncread(file,'v10'),2) , [2 1 3]);

height_mean = sum(tim.Hub_height.*tim.Number_of_turbines) ./ sum(tim.Number_of_turbines);

alpha = (log(w_u100)-log(w_u10)) ./ (log(100)-log(10));
wu_height = w_u100 .* (height_mean/100).^alpha;
alpha = (log(w_v100)-log(w_v10)) ./ (log(100)-log(10));
wv_height = w_v100 .* (height_mean/100).^alpha;

Fu = griddedInterpolant({w_lat,w_lon,datenum(w_time)},wu_height,'linear','nearest');
wu = Fu({g.lat,g.lon,datenum(time)});
Fv = griddedInterpolant({w_lat,w_lon,datenum(w_time)},wv_height,'linear','nearest');
wv = Fv({g.lat,g.lon,datenum(time)});

clear Fu Fv wv_height wu_height w_* alpha

% 
% ws = abs(sqrt(wu.^2+wv.^2));
wdir = angle(wu + 1i*wv);


%% ratio
g_ratio = predict(compactmdl, [repmat(g.lat2D(g.latlonmask),g.nt,1) repmat(g.lon2D(g.latlonmask),g.nt,1) g.NNT(:) repelem(datenum(g.time-g.time(1)),g.nlm,1) ones(numel(g.NNT(:)),2)]);
g_ratio = reshape(g_ratio,size(g.NNT));
g_ratioll = nan(g.nlat,g.nlon,g.nt);
g_ratioll(repmat(g.latlonmask,1,1,g.nt)) = g_ratio;

%% bird density template
template15m = nan(g.nlat,g.nlon,g.nt);
template15m(repmat(g.latlonmask,1,1,g.nt)) = 1e-6 .* sqrt(guv.u_est.^2 + guv.v_est.^2)*60*60;

template1h = turbineSweptAreaMap .* abs(cos(wdir-birdDir));



%% Loop

for ii=0:5
    load(['../Density_simulationMap_reassemble/Density_simulationMap_reassemble_' num2str(ii)]);
    for i_real = 1:size(real_dens,4)
        % rar(:, (ii-1)*100+i_real) = reshape(nansum(nansum(birdAtRisk_t .* real_dens(:,:,:,i_real),1),2),1,[]);
        
        dens = real_dens(:,:,:,i_real);
        % dens = nan(g.nlat,g.nlon,g.nt);dens(repmat(g.latlonmask,1,1,g.nt)) = gd.dens_est;
        
        % Compute ratio uncertainty
        g_ratio_std = polyval(p_abserror,dens*1000);
        ratio_sim = g_ratioll+normrnd(0,g_ratio_std);
        
        % symetry of error on the positive side to maintain not below zero
        ratio_sim(ratio_sim<0)=0-ratio_sim(ratio_sim<0);

        % Compute birdflow
        birdFlow_15min = nan(g.nlat,g.nlon,g.nt);
        birdFlow_15min(repmat(g.latlonmask,1,1,g.nt)) = dens .* ratio_sim .* template15min; 

        birdFlow = nan(g.nlat,g.nlon,numel(time));
        for i=1:numel(time)
            birdFlow(:,:,i) = nanmean(birdFlow_15min(:,:,Locb==i),3);
        end

        birdAtRisk = birdFlow .* template1hr;

    end
end



