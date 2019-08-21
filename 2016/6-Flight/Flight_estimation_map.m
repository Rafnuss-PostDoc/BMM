%% Map Estimation

%% 1. Load data
%

clear all; 
load('../1-Cleaning/data/dc_corr.mat'); 
load('./data/Flight_inference.mat'); 
load('../4-Estimation/data/Density_estimationMap','g');
addpath('./functions/'); 

ndc = numel(dc);
dclat = [dc.lat]';
dclon = [dc.lon]';

clear dc

%% Matrix of Distance
Dtime = squareform(pdist(data.time));
Ddist_sm = squareform(pdist([dclat dclon], @lldistkm));
Ddist = Ddist_sm(data.radar,data.radar);

%% Kriging U
tic; Caa = uv.cov.Cu(Ddist,Dtime);
Dtime_sm = single(pdist2(data.time, datenum(g.time')));
Ddist_sm = single(pdist2([dclat dclon], [g.lat2D(:) g.lon2D(:)], @lldistkm));
guv_ut_est = nan(size(g.mask_day));
guv_ut_sig = nan(size(g.mask_day));
for i_t=1:g.nat
    disp(i_t); toc
    mask_r = data.day==i_t;
    mask_g = g.mask_day & repmat(reshape(g.dateradar==i_t,1,1,[]),g.nlat,g.nlon,1) & repmat(g.latlonmask,1,1,g.nt);
    [tmp_i,tmp_j,tmp_k] = ind2sub(size(g.time3D),find(mask_g));
    
    Dtime_sf = Dtime_sm(mask_r,tmp_k);
    Ddist_sf = Ddist_sm(data.radar(mask_r), sub2ind(size(g.lat2D),tmp_i,tmp_j));
    
    Cab = uv.cov.Cu(Ddist_sf,Dtime_sf);

    Lambda = inv(Caa(mask_r,mask_r))*Cab;

    guv_ut_est(mask_g) = Lambda' * uv.ut(mask_r);
    guv_ut_sig(mask_g) = sqrt(sum(uv.cov.parm_u(1:2)) - sum(Lambda.*Cab)); 
    
end

%% Kriging V
tic; Caa = uv.cov.Cv(Ddist,Dtime);
Dtime_sm = single(pdist2(data.time, datenum(g.time')));
Ddist_sm = single(pdist2([dclat dclon], [g.lat2D(:) g.lon2D(:)], @lldistkm));
guv_vt_est = nan(size(g.mask_day));
guv_vt_sig = nan(size(g.mask_day));
for i_t=1:g.nat
    disp(i_t); toc
    mask_r = data.day==i_t;
    mask_g = g.mask_day & repmat(reshape(g.dateradar==i_t,1,1,[]),g.nlat,g.nlon,1) & repmat(g.latlonmask,1,1,g.nt);
    [tmp_i,tmp_j,tmp_k] = ind2sub(size(g.time3D),find(mask_g));
    
    Dtime_sf = Dtime_sm(mask_r,tmp_k);
    Ddist_sf = Ddist_sm(data.radar(mask_r), sub2ind(size(g.lat2D),tmp_i,tmp_j));
    
    Cab = uv.cov.Cv(Ddist_sf,Dtime_sf);

    Lambda = inv(Caa(mask_r,mask_r))*Cab;

    guv_vt_est(mask_g) = Lambda' * uv.vt(mask_r);
    guv_vt_sig(mask_g) = sqrt(sum(uv.cov.parm_v(1:2)) - sum(Lambda.*Cab)); 
    
end

%% Back-transform
%

guv.u_est = guv_ut_est *uv.t.std(1) + uv.t.mean(1);
guv.v_est = guv_vt_est *uv.t.std(2) + uv.t.mean(2);

guv.u_sig = sqrt(guv_ut_sig.^2 * uv.t.std(1)^2 );
guv.v_sig = sqrt(guv_vt_sig.^2 * uv.t.std(2)^2 );


%% Figure
u=guv.u_est;
v=guv.v_est;
u_isnan=single(~isnan(guv.u_est));
v_isnan=single(~isnan(guv.v_est));
u(isnan(u))=0;
v(isnan(v))=0;

rzd=1/4;
lat2D_res=imresize(g.lat2D,rzd);
lon2D_res=imresize(g.lon2D,rzd);

h=figure(2);  
worldmap([min(g.lat) max(g.lat)], [min(g.lon) max(g.lon)]);  
filename='data/Flight_estimationMap';
mask_fullday = find(~reshape(all(all(isnan(guv.u_est),1),2),g.nt,1));
Frame(numel(mask_fullday)-1) = struct('cdata',[],'colormap',[]); 
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow('worldrivers.shp','Color', 'blue')
set(gcf,'color','w');
for i_t = 1:numel(mask_fullday)

    u_res = imresize(u(:,:,mask_fullday(i_t)),rzd);
    u_isnan_res = imresize(u_isnan(:,:,mask_fullday(i_t)),rzd);
    u_res=u_res./u_isnan_res;
    u_res(u_isnan_res<0.5)=nan;

    v_res = imresize(v(:,:,mask_fullday(i_t)),rzd);
    v_isnan_res = imresize(v_isnan(:,:,mask_fullday(i_t)),rzd);
    v_res=v_res./v_isnan_res;
    v_res(v_isnan_res<0.5)=nan;
    
    hsurf=quiverm(lat2D_res,lon2D_res,u_res,v_res,'k');
    
    drawnow
    title(datestr(g.time(mask_fullday(i_t)))); drawnow;
    Frame(i_t) = getframe(h);
    [imind,cm] = rgb2ind(frame2im(Frame(i_t)),256); 
    if i_t == 1
        imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append','DelayTime',0.1);
    end
    delete(hsurf);
end

v=VideoWriter([filename '.avi']);
v.FrameRate = 4;
v.Quality = 75;
open(v);
writeVideo(v, Frame);
close(v);

%%  Save
%
save('data/Flight_estimationMap','g','guv','-v7.3')
% load('data/Flight_estimationMap')