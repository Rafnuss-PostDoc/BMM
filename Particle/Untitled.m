
load('../2018/data/Density_estimationMap','g','gd');
load('../2018/data/Flight_estimationMap','guv');
load('../2018/data/SinkSource.mat')
addpath('../2018/functions/')

vy = nan(g.nlat, g.nlon, gext.nt); vx=vy; rho = vy;
tmp = repmat(g.latlonmask,1,1,gext.nt); tmp(:,:,ismember(gext.time,g.day-.5))=false;
vy(tmp) = guv.v_est/1000*60*60; % m/s -> km/h (+) north, (-) south
vx(tmp) = guv.u_est/1000*60*60; 

tmp2 = gd.dens_est; % bird/km^2 .* repmat(area(g.latlonmask),1,g.nt); % bird
tmp2(g.rain(repmat(g.latlonmask,1,1,g.nt))>data.mask_rain_thr) = 0;
rho(tmp) = tmp2;

dt = 15/60; 
dy = lldistkm([g.lat(1) g.lon(1)],[g.lat(2) g.lon(1)]);
dx = lldistkm([g.lat2D(:,1) g.lon2D(:,1)],[g.lat2D(:,1) g.lon2D(:,2)]);

%%
F_vlat=griddedInterpolant({1:gext.nlat,1:gext.nlon,1:gext.nt},padarray(vy,[1 1 0],nan),'nearest');
F_vlon=griddedInterpolant({1:gext.nlat,1:gext.nlon,1:gext.nt},padarray(vx,[1 1 0],nan),'nearest');
F_land=griddedInterpolant({1:gext.nlat,1:gext.nlon,1:gext.nt-1},padarray(-W.landing./rho(:,:,1:end-1),[1 1 0],1),'nearest');
F_out = griddedInterpolant({1:gext.nlat,1:gext.nlon},padarray(double(g.latlonmask),[1 1],0),'nearest','none');
% Number of birds to simulate: 1 for every 1000.
% nb_bird=1000*ones(g.nat-1,1); % Fixte number
nb_bird=round((Ts.Fout.day.entering + Ts.day.takingoff));



for i_day=1:g.nat
    % Find the time of the night
    idt=find(gext.day_id(1:end-1)==i_day);
    
    tmp = Fout.entering(:,:,idt);
    tmp(2:end-1,2:end-1,:) = tmp(2:end-1,2:end-1,:)+W.takingoff(:,:,idt);
    
    y = randsample(numel(tmp), nb_bird(i_day), true, tmp(:));
    
    y_lat = nan(nb_bird(i_day), numel(idt)); 
    y_lon = y_lat;
    [y_lat(:,2), y_lon(:,2), y_t] = ind2sub(size(tmp),y);
    
    y_t=y_t+1;
    
    for i_t=2:numel(idt)
        % Find bird in the air
        id = (y_t <= i_t);
     
        % Find if bird lands until next step (no need to move them forward)
        landing = F_land(y_lat(:, i_t), y_lon(:, i_t), repmat(idt(i_t-1),nb_bird(i_day),1))>rand(nb_bird(i_day),1);
        y_t(id & landing)=Inf;
        id(landing)=false;
        
        % Find the speed
        v_lat = F_vlat(y_lat(id, i_t), y_lon(id, i_t), repmat(idt(i_t),sum(id),1)); % km/h
        v_lon = F_vlon(y_lat(id, i_t), y_lon(id, i_t), repmat(idt(i_t),sum(id),1)); % km/h
        
        % Update position
        y_lat(id, i_t+1) = y_lat(id, i_t) + v_lat*dt/dy;
        y_lon(id, i_t+1) = y_lon(id, i_t) + v_lon.*dt./interp1((1:g.nlat)',dx,y_lat(id, i_t),'linear','extrap');
        
        % keep same position for that have not yet left
        y_lat(~id, i_t+1) = y_lat(~id, i_t);
        y_lon(~id, i_t+1) = y_lon(~id, i_t); 

        % Check for errors
        % assert(~any(isnan(y_lon(:, i_t+1))))
        y_lat(isnan(y_lat(:, i_t+1)), i_t+1) = y_lat(isnan(y_lat(:, i_t+1)), i_t);
        y_lon(isnan(y_lon(:, i_t+1)), i_t+1) = y_lon(isnan(y_lon(:, i_t+1)), i_t); 
        
        % Check if bird out
        out = ~(F_out(y_lat(:, i_t+1), y_lon(:, i_t+1))==1);
        y_t(out) = 999;
        
    end
    
    % assert(~any(isnan(y_lon(:))))
    assert(all(~y_t<100))
    % Compute stat
    
    ori = sub2ind([gext.nlon gext.nlat],round(y_lon(:,2)),round(y_lat(:,2)));
    dest = sub2ind([gext.nlon gext.nlat],round(y_lon(:,i_t+1)),round(y_lat(:,i_t+1)));
    [C,~,ic] = unique([ori dest],'rows');
    res{i_day} = [C splitapply(@sum,ones(numel(ic),1),ic)];
end