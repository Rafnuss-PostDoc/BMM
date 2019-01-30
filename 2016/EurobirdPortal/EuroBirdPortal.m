%% Corrected map

% Coordinate
n = [280 226];
dxy = 30000;
x = 0 + dxy*(0:n(1)-1);
y =  720000 + dxy*(n(2)-1:-1:0);
[X,Y]=meshgrid(x,y);

% Convert in lat, lon
fn = 'coord_3035.txt';
fno = 'coord_4326.txt';
dlmwrite(fn,[X(:) Y(:)],'delimiter',' ','precision','%7.f')
dlmwrite(fn,[X(:) Y(:)],'delimiter',' ','precision','%7.f')
status = system(['bash -c "gdaltransform -s_srs EPSG:3035 -t_srs EPSG:4326 <  ' fn '  >  ' fno '"']);
LL = dlmread(fno);
LAT = reshape(LL(:,2),n(2),n(1));
LON = reshape(LL(:,1),n(2),n(1));

% laod data of corrected map
da=38:41;
OENOEN = nan(n(2),n(1),numel(da));
TURPHI = nan(n(2),n(1),numel(da));
for i=1:numel(da)
    OENOEN(:,:,i)=reshape(importdata(['data\EBP_files\maps\OENOEN_2016_' num2str(da(i)) '.asc']),n(1),n(2))';
    TURPHI(:,:,i)=reshape(importdata(['data\EBP_files\maps\TURPHI_2016_' num2str(da(i)) '.asc']),n(1),n(2))';
end
OENOEN(OENOEN==-9999)=nan;
TURPHI(TURPHI==-9999)=nan;



% Load data from BMM
load('data/Density_estimationMap','g')
c = (g.dens_est); % bird/km^2
load('data/FlightSpeed_estimationMap','g')
speed = g.dens_est * 60*60/1000; % convert to km/hr

% MTR per night
MTR = c .* speed;
MTR_day=nan(g.nlat,g.nlon,g.nat);
for i_t=1:g.nat-1
    idt=g.atime(i_t)-0.5<datenum(g.time) & g.atime(i_t)+0.5>datenum(g.time);
    MTR_day(:,:,i_t)=nanmean(MTR(:,:,idt),3); 
end
MTR_day(~repmat(g.latlonmask,1,1,g.nat))=nan;

% load coastlines;
% figure('position',[0 0 1000 600]); 
% for i_t = 1:g.nat-1
%     subplot(4,6,i_t); hold on; worldmap([g.lat(1) g.lat(end)], [g.lon(1) g.lon(end)]);
%     surfm(g.lat2D,g.lon2D,MTR_day(:,:,i_t)); caxis([0 3000])
%     plotm(coastlat, coastlon,'k')
% end

% Average of the period of time
dds=datetime({'16-Sep-2016 12:00','23-Sep-2016 12:00','30-Sep-2016 12:00','7-Oct-2016 12:00','23-Oct-2016 12:00'},'InputFormat','dd-MMM-yyy HH:mm');
MTR_week=nan(g.nlat,g.nlon,numel(dds)-1);
for i_t=1:numel(dds)-1
    idt=datenum(dds(i_t))<g.atime & datenum(dds(i_t+1))>g.atime;
    MTR_week(:,:,i_t)=nanmean(MTR_day(:,:,idt),3); 
end




% Plot 
load coastlines;
figure(4);clf; colormap('winter')
for i_t = 1:numel(da)
    %subplot(3,numel(da),i_t); 
    worldmap([43 60], [-5 20]);gridm('off');
    geoshow('landareas.shp', 'FaceColor', [27 27 33]/255)
    surfm(g.lat2D,g.lon2D,MTR_week(:,:,i_t)); caxis([0 2000])
    plotm(coastlat, coastlon,'k')
    print(['MTR_' num2str(i_t) '.eps'],'-depsc')
     
    %subplot(3,numel(da),i_t+4); 
    worldmap([43 60], [-5 20]);gridm('off');
    geoshow('landareas.shp', 'FaceColor', [27 27 33]/255)
    surfm(LAT,LON,OENOEN(:,:,i_t)); caxis([0 1])
    plotm(coastlat, coastlon,'k') 
    print(['OENOEN_' num2str(i_t) '.eps'],'-depsc')
    
    %subplot(3,numel(da),i_t+8);
    worldmap([43 60], [-5 20]); gridm('off');
    geoshow('landareas.shp', 'FaceColor', [27 27 33]/255)
    surfm(LAT,LON,TURPHI(:,:,i_t)); caxis([0 1])
    plotm(coastlat, coastlon,'k') 
    set(gca,'YtickLabel',[]);
    print(['TURPHI_' num2str(i_t) '.eps'],'-depsc')
end



%% Load data frp, EBP 
% NO USED: corrected map used instead
R = readtable('data\EBP_files\Records_2016.csv');
E = readtable('data\EBP_files\Events_2016.csv');
load('data/RE')
load coastlines.mat;



% Record and Event restructure
RE.xy = unique([R.MAPX R.MAPY],'rows');
RE.t = unique(R.DATE_WEEK);
RE.XY = repmat(RE.xy,1,1,numel(RE.t));
RE.T = repmat(reshape(RE.t,1,1,[]),size(RE.xy,1),size(RE.xy,2),1);

% Create sctructure from variable of table and fill them
Evn = E.Properties.VariableNames(7:end);
for i_vn=1:numel(Evn)
    E.(char(Evn(i_vn)))(isnan(E.(char(Evn(i_vn)))))=0; % replace nan by 0
    RE.(string(Evn(i_vn)))= zeros(size(RE.xy,1),size(RE.t,1));
end
for e=1:size(E,1)
    id_xy = find(all(RE.xy==[E.MAPX(e) E.MAPY(e)],2));
    id_t = find(RE.t==E.DATE_WEEK(e));
    for i_vn=1:numel(Evn)
        RE.(string(Evn(i_vn)))(id_xy,id_t) = RE.(string(Evn(i_vn)))(id_xy,id_t)+E.(char(Evn(i_vn)))(e);
    end
end

RE.species = sort(unique(R.SPECIES_CD));
Rvn = R.Properties.VariableNames(8:end);
for i_sp=1:numel(RE.species)
    R.(char(Rvn(i_vn)))(isnan(R.(char(Rvn(i_vn)))))=0; % replace nan by 0
    RE.(string(RE.species(i_sp)))=struct();
    for i_vn=1:numel(Rvn)
        RE.(string(RE.species(i_sp))).(string(Rvn(i_vn))) = zeros(size(RE.xy,1),size(RE.t,1));
    end
end
for r=1:size(R,1)
    id_xy = find(all(RE.xy==[R.MAPX(r) R.MAPY(r)],2));
    id_t = find(RE.t==R.DATE_WEEK(r));
    for i_vn=1:numel(Rvn)
        RE.(string(R.SPECIES_CD(r))).(string(Rvn(i_vn)))(id_xy,id_t) = RE.(string(RE.species(i_sp))).(string(Rvn(i_vn)))(id_xy,id_t) + R.(char(Rvn(i_vn)))(r);
    end
end


% Convert coordinate
fn = 'coord_3035.txt';
fno = 'coord_4326.txt';
dlmwrite(fn,RE.XY(:,:,1),'delimiter',' ','precision','%7.f')
status = system(['bash -c "gdaltransform -s_srs EPSG:3035 -t_srs EPSG:4326 <  ' fn '  >  ' fno '"']);
RE.ll = dlmread(fno);
RE.ll = RE.ll(:,1:2);
RE.LL = repmat(RE.ll,1,1,numel(RE.t));

figure; worldmap([min(RE.ll(:,2)) max(RE.ll(:,2))],[  min(RE.ll(:,1)) max(RE.ll(:,1))]);
scatterm(RE.ll(:,2),RE.ll(:,1),'.k')

% save('data/RE','RE')


% Explore data
figure(3);
for i=1:numel(Evn)
    subplot(3,2,i); histogram(RE.(string(Evn{i})));
    %set(gca,'xscale','log'); 
    xlabel(string(Evn{i}))
end


figure(4);
for i=1:numel(Evn)
    subplot(3,2,i); histogram(RE.OENOEN.(string(Rvn{i})));
    %set(gca,'xscale','log'); 
    xlabel(string(Evn{i}))
end

figure; subplot(1,2,1); hold on;
worldmap([min(RE.ll(:,2)) max(RE.ll(:,2))],[  min(RE.ll(:,1)) max(RE.ll(:,1))]); plotm(coastlat, coastlon,'k')
scatterm(RE.ll(:,2),RE.ll(:,1),[],nansum(RE.COMP_LISTS,2),'filled')
colorbar; title('Number of Complete List')
subplot(1,2,2); hold on;
worldmap([min(RE.ll(:,2)) max(RE.ll(:,2))],[  min(RE.ll(:,1)) max(RE.ll(:,1))]); plotm(coastlat, coastlon,'k')
scatterm(RE.ll(:,2),RE.ll(:,1),[],nansum(RE.CASUAL_EV,2),'filled')
colorbar; title('Number of Casual Report')

figure; subplot(1,2,1); hold on;
worldmap([min(RE.ll(:,2)) max(RE.ll(:,2))],[  min(RE.ll(:,1)) max(RE.ll(:,1))]); plotm(coastlat, coastlon,'k')
scatterm(RE.ll(:,2),RE.ll(:,1),[],nansum(RE.OENOEN.COMP_LISTS,2)./nansum(RE.COMP_LISTS,2),'filled')
colorbar; title('Ratio of complete list with report of wheatear'); caxis([0 1])
subplot(1,2,2); hold on;
worldmap([min(RE.ll(:,2)) max(RE.ll(:,2))],[  min(RE.ll(:,1)) max(RE.ll(:,1))]); plotm(coastlat, coastlon,'k')
scatterm(RE.ll(:,2),RE.ll(:,1),[],nansum(RE.OENOEN.CASUAL_EV,2)./nansum(RE.CASUAL_EV,2),'filled')
colorbar; title('Ratio of casual event with report of wheatear'); caxis([0 1])

figure(3); colormap('winter')
for i=1:4
    ax =subplot(1,4,i); hold on; 
    % worldmap([min(RE.ll(:,2)) max(RE.ll(:,2))],[  min(RE.ll(:,1)) max(RE.ll(:,1))]); 
    worldmap([41 60], [-5 20]);
    geoshow(ax, 'landareas.shp', 'FaceColor', [0.5 0.7 0.5])
    geoshow('worldrivers.shp','Color', 'blue')
    pst =(RE.OENOEN.MAX_COUNT(:,i)>0);
    scatterm(RE.ll(~pst,2),RE.ll(~pst,1),4,[0.4 0.4 0.4],'filled');
    pst=find(pst);
    [a,id] = sort(RE.OENOEN.MAX_COUNT(pst,i));
    sz = log(a)*2+1;
    scatterm(RE.ll(pst(id),2),RE.ll(pst(id),1),sz+5,(a),'filled')
    caxis([5 200])
    % plotm(coastlat, coastlon,'k')
end






