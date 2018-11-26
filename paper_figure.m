%% Load data
clear all; load('./data/dc_corr.mat'); load coastlines; addpath('functions/'); warning('off'); load('data/Density_modelInf.mat'); load('data/Density_estimationMap')


%% Figure 2: Radar availability
load('./data/d_all.mat')

for i_d=1:numel(d)
    d_i_d = datenum(getabstime(d(i_d).dens));
    n(i_d) = sum(unique(round(d_i_d))>datenum(start_date) & unique(round(d_i_d))<datenum(end_date));
end

figure('position',[0 0 800 600]);
worldmap([min([d.lat]) max([d.lat])], [min([d.lon]) max([d.lon])]);  
files={'gt30w020n90','gt30e020n40','gt30w020n40','gt30e020n90'};
for i=1:numel(files)
    [X,cmap,R] = geotiffread(['C:\Users\rnussba1\Documents\MATLAB\' files{i}]);
    geoshow(double(X),cmap,'DisplayType','texturemap')
end
demcmap(X)
% export_fig 'figure/paper/radarsNetwork_1.eps' -eps

figure('position',[0 0 800 600]);
worldmap([min([d.lat]) max([d.lat])], [min([d.lon]) max([d.lon])]);  
% plotm(coastlat, coastlon,'k','linewidth',2)
h1=scatterm([d.lat],[d.lon],100,n,'filled','MarkerEdgeColor','k'); 
h2=scatterm([dc.lat],[dc.lon],100,'MarkerEdgeColor','r');
c=colorbar('southoutside');c.Label.String='Number of day with data between 19-sept. to 10-oct.';
%textm([d.lat]-.4,[d.lon],{d.name},'HorizontalAlignment','center')
legend([h1,h2],{'Radars of network','Radars used in the study'})

% export_fig 'figure/paper/radarsNetwork_2.eps' -eps


%% Stats Radars
load('./data/d.mat')
ddens=[];
for i=1:numel(d)
    ddens=[ddens;d(i).dens(:)];
end
figure; hold on; 
histogram(ddens,'EdgeColor','none','FaceAlpha',1);
histogram(data.dens,'EdgeColor','none','FaceAlpha',1);
set(gca, 'YScale', 'log'); axis([0 600 1 10^5]); %set(gca, 'XScale', 'log')
% export_fig 'figure/paper/hist_raw_data.eps' -eps

S=nan(numel(dc),1);
for i=1:numel(dc)
    S(i)=numel(data.dens(data.i_r==i));
end
figure; hold on; 
histogram(S,40,'EdgeColor','none','FaceAlpha',1);

Ddist_sf=squareform(pdist([[dc.lat]' [dc.lon]'],@lldistkm));
Ddist_sf(Ddist_sf==0)=Inf;
figure; hold on; 
histogram(min(Ddist_sf),'EdgeColor','none','FaceAlpha',1);
plot(min(Ddist_sf),0,'.k')

%% Figure 1: Illustration of a radar

fig0=figure('position',[0 0 1000 400]); clf; set(gcf, 'Color', 'w');
i_d=find(strcmp({dc.name},'bezav'));

subplot(4,1,1); hold on;
imagesc(datenum(d(i_d).time), d(i_d).interval*(1/2:double(d(i_d).levels)), d(i_d).DBZH','AlphaData',~isnan(d(i_d).DBZH'))
xlabel('Date'); ylabel('Altitude [m]'); c=colorbar; c.Label.String='Reflectivity';
datetick('x'); axis([datenum(start_date) datenum(end_date-1) 0 5000]); set(gca, 'YDir', 'normal');ylim([0 4000])

subplot(4,1,2); hold on;
imagesc(datenum(d(i_d).time), d(i_d).interval*(1/2:double(d(i_d).levels)), log(d(i_d).dens)','AlphaData',~isnan(d(i_d).dens'))
plot([datenum(start_date) datenum(end_date)],[d(i_d).height d(i_d).height],'r','linewidth',2); caxis([-5 5])
xlabel('Date'); ylabel('Altitude [m]'); c=colorbar; c.Label.String='Bird Density [bird/m^3]';
datetick('x'); axis([datenum(start_date) datenum(end_date-1) 0 5000]); set(gca, 'YDir', 'normal');ylim([0 4000])

subplot(4,1,3); hold on;
imagesc(datenum(d(i_d).time), d(i_d).interval*(1/2:double(d(i_d).levels)), log(dc(i_d).dens)','AlphaData',~isnan(d(i_d).dens'))
plot([datenum(start_date) datenum(end_date)],[d(i_d).height d(i_d).height],'r','linewidth',2); caxis([-5 5])
xlabel('Date'); ylabel('Altitude [m]'); c=colorbar; c.Label.String='Bird Density [bird/m^3]';
datetick('x'); axis([datenum(start_date) datenum(end_date-1) 0 5000]); set(gca, 'YDir', 'normal');ylim([0 4000])

subplot(4,1,4); hold on;
plot(data.time(data.i_r==i_d), data.dens(data.i_r==i_d),'.k'); axis tight
% export_fig 'figure/paper/Illustration_radar_DHZ_DENS3.eps' -eps


%% Figure 3: Features illustration

fig1=figure('position',[0 0 1200 500]);clf;
ncol=6; nrow=3;
ax=subplot(nrow,ncol,[1 ncol*(nrow-1)+2]); hold on; 
worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]);  
geoshow(ax, 'landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow('worldrivers.shp','Color', 'blue')
%plotm(coastlat, coastlon,'k')
S=zeros(numel(dc), 1);
for i=1:numel(dc)
    S(i)=mean(data.dens(data.i_r==i));
end
scatterm([dc.lat],[dc.lon],100,S,'filled','MarkerEdgeColor','k'); colorbar('southoutside');

co=get(gca,'ColorOrder');u=0;
co=[51,34,136;136,204,238;68,170,153;17,119,51;153,153,51;221,204,119;204,102,119;136,34,85;170,68,153]/255;
radars={{'fikuo','fiuta','fivim'},{'denhb','bewid','frave'},{'frgre','frbor','frmom'}}; %
for i=1:numel(radars)
    for j=1:numel(radars{i})
        u=u+1;
        i_d=find(strcmp({dc.name},radars{i}{j}));
        subplot(nrow,ncol,ncol*i+[-3 -1]); hold on; box on;
        if j==1
            X=[unique(dc(i_d).sunrise)' unique(dc(i_d).sunset)'];
            for s=1:size(X,1)-1
                fill([X(s,1) X(s,1) X(s+1,2) X(s+1,2)],[0 200 200 0 ],[252, 255, 196]./255,'EdgeColor','none')
            end
        end
        plot(datetime(data.time(data.i_r==i_d),'ConvertFrom','datenum'), data.dens(data.i_r==i_d),'.','Color',co(u,:));
        ylabel('Bird density [bird/m^2]')
        xlim([start_date end_date]); ylim([0 200])
        subplot(nrow,ncol,ncol*i); hold on; box on
        plot(datetime(data.time(data.i_r==i_d),'ConvertFrom','datenum'), data.dens(data.i_r==i_d),'.','Color',co(u,:));
        xlim(datetime({'26/09/2016 12:00' '28/09/2016 12:00'}) );
        if j==numel(radars{i})
            yl=get(gca,'ylim');
            for s=1:size(X,1)-1
                fill([X(s,1) X(s,1) X(s+1,2) X(s+1,2)],[yl fliplr(yl)],[252, 255, 196]./255,'EdgeColor','none')
            end
        end
        scatterm(ax,dc(i_d).lat,dc(i_d).lon,100,'MarkerEdgeColor',co(u,:),'linewidth',3);
    end
end

% export_fig 'figure/paper/conceptual_model_3each.eps' -eps

%% Figure 4: Mathematical model
fig2=figure('position',[0 0 1000 400]);clf; clf; hold on; box on;
i_d=find(strcmp({dc.name},'frmom'));
%co=get(gca,'ColorOrder');

t=39;
A=[53 36 42];
d=[594 663 732];

plot(datetime({'26/09/2016 12:00' '29/09/2016 12:00'}),[t t],'.-k');
poly = [1^2 1 1;44^2 44 1;22^2 22 1] \ [-30;-30;10];

for i=1:numel(A)
    d_d=data.dens(data.dateradar==d(i));
    d_a=datetime(data.time(data.dateradar==d(i)),'ConvertFrom','datenum');
    plot([d_a(1) d_a(end)],[A(i) A(i)],'--k');
    plot(linspace(d_a(1),d_a(end),44),A(i)+[(1:44)'.^2 (1:44)' ones(44,1)]*poly,'k')
    % fill([linspace(d_a(1),d_a(end),44) d_a(end) d_a(1)],[ A(i)+[(1:44)'.^2 (1:44)' ones(44,1)]*poly ; t ;t],[.95 .95 .95],'edgecolor','none')
    c_a=interp1(linspace(d_a(1),d_a(end),44),A(i)+[(1:44)'.^2 (1:44)' ones(44,1)]*poly,d_a);
    for ii=1:numel(c_a)
        plot([d_a(ii) d_a(ii)],[c_a(ii) d_d(ii)],'Color',[.6 .6 .6])
    end
end
plot(datetime(data.time(data.i_r==i_d),'ConvertFrom','datenum'), data.dens(data.i_r==i_d),'.k')
xlim(datetime({'26/09/2016 12:00' '29/09/2016 12:00'}) ); ylabel('Bird density [bird/m^2]');

% export_fig 'figure/paper/mathematical_model.eps' -eps


%% Figure 5: Kriging estimation
latlonmask = find(g.latlonmask);
Gneiting = @(dist,time,range_dist,range_time,delta,gamma,beta) 1./( (time./range_time).^(2.*delta) +1 ) .* exp(-(dist./range_dist).^(2.*gamma)./((time./range_time).^(2.*delta) +1).^(beta.*gamma) );
wradius=4;

i_l=g.nlm/3;
id1=find(gampli.Ddist_sf(:,latlonmask(i_l))<ampli.cov.parm(3)*wradius);
i_t=g.nat/2;
id2= id1(gampli.Dtime_sf(id1,i_t)<ampli.cov.parm(4)*wradius);
Cab = ampli.cov.parm(2).*Gneiting(gampli.Ddist_sf(id2,latlonmask(i_l)), gampli.Dtime_sf(id2,i_t), ampli.cov.parm(3), ampli.cov.parm(4), ampli.cov.parm(5), ampli.cov.parm(6), ampli.cov.parm(7));
lambda = ampli.cov.C(id2,id2)  \  Cab;
An_est = lambda' * ampli.An(id2);
An_sig = sqrt(ampli.cov.parm(1) + ampli.cov.parm(2) - lambda' * Cab);

figure; hold on;
idisnan=find(~ampli.isnanA);
[I_d, I_t]=ind2sub(size(ampli.A),idisnan(id2));
plot3( coastlon,coastlat, min(g.atime(I_t))*ones(size(coastlon)),'k','linewidth',1)
fill3( coastlon,coastlat, min(g.atime(I_t))*ones(size(coastlon)),'k')
plot3( coastlon,coastlat, 0.5*(max(g.atime(I_t))+min(g.atime(I_t)))*ones(size(coastlon)),'k','linewidth',0.5)
plot3( coastlon,coastlat, max(g.atime(I_t))*ones(size(coastlon)),'k','linewidth',.5)
dc_lat=[dc.lat]; dc_lon=[dc.lon];
scatter3(g.lon2D(latlonmask(i_l)),g.lat2D(latlonmask(i_l)),g.atime(i_t),'r','filled')
scatter3(dc_lon(I_d),dc_lat(I_d),g.atime(I_t),[],log(lambda),'filled')
ylim([min(dc_lat(I_d))-3 max(dc_lat(I_d))])
xlim([min(dc_lon(I_d))-3 max(dc_lon(I_d))+2]); 
zlim([min(g.atime(I_t)) max(g.atime(I_t))])
view(3); box on; ylabel('Latitude'); xlabel('Longitude'); zlabel('Time '); datetick('z','dd-mmm','keeplimits')
colorbar('southoutside');
%set(gca,'BoxStyle','full'); 
grid on
% export_fig 'figure/paper/kriging.eps' -eps



%% Figure 6: Result Inferance

fig6=figure('position',[0 0 1400 700]); 
subplot(2,4,[1 2]); hold on;
worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]); 
[LAT,LON] = meshgrid(min([dc.lat]):max([dc.lat]),min([dc.lon]):max([dc.lon]));
surfm(LAT,LON,LAT*trend.p(1)+LON*trend.p(2)+trend.p(3))
plotm(coastlat, coastlon,'k');
geoshow('worldrivers.shp','Color', 'blue')
colorbar;
%print -depsc2 figure/paper/result_calibration_1.eps;

%fig6=figure('position',[0 0 1400 700]); 
subplot(2,4,[3 4]); hold on;
x=(-1:.01:1)';
y=polyval( curve.p, x);
z=sqrt(polyval( res.p, x));
fill([x ; flipud(x)], [y-3*z ; flipud(y+3*z)],[.9 .9 .9],'EdgeColor','none')
fill([x ; flipud(x)], [y-2*z ; flipud(y+2*z)],[.8 .8 .8],'EdgeColor','none')
fill([x ; flipud(x)], [y-1*z ; flipud(y+1*z)],[.7 .7 .7],'EdgeColor','none')
plot(x,y,'-k','linewidth',2)
xlabel('Normalized Night Time (NNT)'); ylabel('Normalized Bird density');  box on; %set(gca,'YScale','log')
%print -depsc2 figure/paper/result_calibration_2.eps;

%fig6=figure('position',[0 0 1400 700]); 
Gneiting = @(dist,time,range_dist,range_time,delta,gamma,beta) 1./( (time./range_time).^(2.*delta) +1 ) .* exp(-(dist./range_dist).^(2.*gamma)./((time./range_time).^(2.*delta) +1).^(beta.*gamma) );
tmpD = ampli.cov.d(1):1:ampli.cov.d(end);
tmpT = ampli.cov.t(1):0.1:ampli.cov.t(end);
subplot(2,4,5);
tmpCOV=ampli.cov.parm(2).*Gneiting(tmpD, zeros(size(tmpD)), ampli.cov.parm(3), ampli.cov.parm(4), ampli.cov.parm(5), ampli.cov.parm(6), ampli.cov.parm(7));
tmpCOV(1)=tmpCOV(1)+ampli.cov.parm(1);
plot(tmpD,tmpCOV,'-k','linewidth',2);
xlabel('Distance [km]'); ylabel('Covariance'); box on; xlim([0 2000])
subplot(2,4,6);
tmpCOV=ampli.cov.parm(2).*Gneiting(zeros(size(tmpT)), tmpT, ampli.cov.parm(3), ampli.cov.parm(4), ampli.cov.parm(5), ampli.cov.parm(6), ampli.cov.parm(7));
tmpCOV(1)=tmpCOV(1)+ampli.cov.parm(1);
plot(tmpT,tmpCOV,'-k','linewidth',2);
xlabel('Time [Days]'); ylabel('Covariance'); box on;  xlim([0 10])

%fig6=figure('position',[0 0 1400 700]); 
tmpD = res.cov.d(1):1:res.cov.d(end);
tmpT = res.cov.t(1):0.01:res.cov.t(end);
subplot(2,4,7);
tmpCOV=res.cov.parm(2).*Gneiting(tmpD, zeros(size(tmpD)), res.cov.parm(3), res.cov.parm(4), res.cov.parm(5), res.cov.parm(6), res.cov.parm(7));
tmpCOV(1)=tmpCOV(1)+res.cov.parm(1);
plot(tmpD,tmpCOV,'-k','linewidth',2);
xlabel('Distance [km]'); ylabel('Covariance'); box on; xlim([0 1000])
subplot(2,4,8);
tmpCOV=res.cov.parm(2).*Gneiting(zeros(size(tmpT)), tmpT, res.cov.parm(3), res.cov.parm(4), res.cov.parm(5), res.cov.parm(6), res.cov.parm(7));
tmpCOV(1)=tmpCOV(1)+res.cov.parm(1);
plot(tmpT,tmpCOV,'-k','linewidth',2);
xlabel('Time [Days]'); ylabel('Covariance'); box on;  xlim([0 0.4])
%print -depsc2 figure/paper/result_calibration_4.eps;

%% Figure 7: Estimation Map grid

fig7=figure('position',[0 0 1400 700]); 
worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]); 
a=ones(size(g.lat2D));
a(g.latlonmask)=0;
b=pcolorm(g.lat2D,g.lon2D,a,'FaceColor',[0 0 51]./255,'FaceAlpha',.3);
plotm(coastlat, coastlon,'k');
geoshow('worldrivers.shp','Color', 'blue')
colorbar;


%% Figure 7: Result Estimation Map

i_t=1+4*12+4*24*12 + [-4*5:4:4*5];
i_lat=[5 45 120];
i_lon=[25 80 165];
ymax=[500 300 100];

c_axis=[-1 2];


figure('position',[0 0 1400 1400]); hold on;
plot3( coastlon,coastlat, datenum(g.time(1))*ones(size(coastlon)),'k','linewidth',2)
%h=slice(g.lon3D,g.lat3D,datenum(g.time3D),log10(g.dens_est),[],[],datenum(g.time(1+4*12:4*24:end)));
%set(h,'edgecolor','none'); grid on;
scatter3(data.lon,data.lat,datenum(data.time),[],data.denstrans,'filled')

ylim([g.lat(1) g.lat(end)])
xlim([g.lon(1) g.lon(end)]); 
zlim(datenum([g.time(1) g.time(end)]));
view(3); box on; ylabel('Latitude'); xlabel('Longitude'); zlabel('Time '); datetick('z','dd-mmm','keeplimits')

% for it=1:numel(i_t)
%     plot3([g.lon(1) g.lon(1) g.lon(end) g.lon(end) g.lon(1)],[g.lat(1) g.lat(end) g.lat(end) g.lat(1) g.lat(1)], datenum(g.time(i_t(it)))*[1 1 1 1 1],'-k','linewidth',2)
% end
% for ill=1:numel(i_lat)
%     plot3( [g.lon(i_lon(ill)) g.lon(i_lon(ill)) ],[g.lat(i_lat(ill)) g.lat(i_lat(ill))], datenum([g.time(1) g.time(end)]),'-r','linewidth',2)
% end
% caxis(c_axis);
%c=colorbar('north'); c.Label.String='Bird Density [bird/m^3]';
% export_fig 'figure/paper/estimation1.eps' -eps


figure('position',[0 0 1400 1400]); hold on;
for it=1:numel(i_t)
    subplot(2,ceil(numel(i_t)/2),it);
    worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]); 
    plotm(coastlat, coastlon,'k');
    geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
    geoshow('worldrivers.shp','Color', 'blue')
    surfm(g.lat2D,g.lon2D,log10(g.dens_est(:,:,i_t(it))));
    for ill=1:numel(i_lat)
        scatterm(g.lat(i_lat(ill)),g.lon(i_lon(ill)),'r','filled')
    end
    %c=colorbar('southoutside'); c.Label.String='Bird Density [bird/m^3]';
    caxis(c_axis);
end

% gtit = datenum(g.time(mask_fullday(i_t)+([-1 0 1])));
% id = find(mean(gtit(1:2))<data.time & mean(gtit(2:3))>data.time);
% if sum(id)>0
%     G = findgroups(data.dateradar(id));
%     hscat=scatterm(splitapply(@mean,data.lat(id),G),splitapply(@mean,data.lon(id),G),[],splitapply(@mean,log(data.dens(id)),G),'filled','MarkerEdgeColor','k');
% end
% for it=1:numel(i_t)
%     %fig7=figure('position',[0 0 1400 1400]); 
%     subplot(nfig,3,(it-1)*6+[3 6]);
%     worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]); 
%     plotm(coastlat, coastlon,'k');
%     geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
%     geoshow('worldrivers.shp','Color', 'blue')
%     surfm(g.lat2D,g.lon2D,log10(g.dens_sig(:,:,i_t(it),5)-g.dens_sig(:,:,i_t(it),3)));
%     for ill=1:numel(i_lat)
%         scatterm(g.lat(i_lat(ill)),g.lon(i_lon(ill)),'r','filled')
%     end
%     %c=colorbar('southoutside'); c.Label.String='Bird Density [bird/m^3]';
%     caxis(c_axis);
%     %export_fig(['figure/paper/estimation' num2str(it) '.eps'],'-eps')
% end


figure('position',[0 0 1400 1400]); hold on; hold on;
co=[51,34,136;136,204,238;68,170,153;17,119,51;153,153,51;221,204,119;204,102,119;136,34,85;170,68,153]/255;
for ill=1:numel(i_lat)
    hold on;
%     a=reshape(g.dens_sig(i_lat(ill),i_lon(ill),:,1),[],1); a(isnan(a))=0;
%     b=reshape(g.dens_sig(i_lat(ill),i_lon(ill),:,7),[],1); b(isnan(b))=0;
%     fill([g.time fliplr(g.time)]', [a ; flipud(b)],[.9 .9 .9],'EdgeColor','none')
%     a=reshape(g.dens_sig(i_lat(ill),i_lon(ill),:,2),[],1); a(isnan(a))=0;
%     b=reshape(g.dens_sig(i_lat(ill),i_lon(ill),:,6),[],1); b(isnan(b))=0;
%     fill([g.time fliplr(g.time)]', [a ; flipud(b)],[.8 .8 .8],'EdgeColor','none')
    a=reshape(g.dens_sig(i_lat(ill),i_lon(ill),:,5),[],1); a(isnan(a))=0;
    b=reshape(g.dens_sig(i_lat(ill),i_lon(ill),:,3),[],1); b(isnan(b))=0;
    fill([g.time fliplr(g.time)]', [a ; flipud(b)],co((ill-1)*3+1,:),'EdgeColor','none','FaceAlpha',.5)
    plot(g.time,reshape(g.dens_sig(i_lat(ill),i_lon(ill),:,4),[],1),'Color',co((ill-1)*3+1,:));
end
for it=1:numel(i_t)
    plot([(g.time(i_t(it))) (g.time(i_t(it)))], [0 ymax(ill)],'k','linewidth',2)
end
xlabel('Date'); ylabel('Bird density'); axis tight;
ylim([0 150]); box on;
%set(gca,'YScale','log')

% export_fig 'figure/paper/estimation2.eps' -eps

%% Figure Integral global
g.dens_est(~g.mask_rain & ~isnan(g.dens_est) )=0;
g_dens_est_avg = nanmean(g.dens_est,3);

g.dens_sig(~repmat(g.mask_rain,1,1,1,nsig) & ~isnan(g.dens_sig) )=0;
g_dens_sig_avg = nanmean(g.dens_sig,3);

figure;
subplot(2,2,1);
worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]); 
plotm(coastlat, coastlon,'k');
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow('worldrivers.shp','Color', 'blue')
surfm(g.lat2D,g.lon2D,log10(g_dens_est_avg));

subplot(2,2,2);
worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]); 
plotm(coastlat, coastlon,'k');
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow('worldrivers.shp','Color', 'blue')
surfm(g.lat2D,g.lon2D,log10(g_dens_sig_avg(:,:,:,5)-g_dens_sig_avg(:,:,:,3)));

subplot(2,2,[3 4]); hold on;
fill([g.time fliplr(g.time)]', [g.sum_sig(:,5) ; flipud(g.sum_sig(:,3))],[.7 .7 .7],'EdgeColor','none','FaceAlpha',.5)
a=g.sum_sig(:,4); a(a==0)=nan;
plot(g.time,a,'k','linewidth',2);
axis tight

%% Figure ... Cross validation

radars={'fiuta','bewid','frgre'}; 
for i_r=1:numel(radars)

    id=data.i_r==find(strcmp({dc.name},radars{i_r}));
    
    subplot(1,numel(radars),i_r);
    
    a=reshape(g.dens_sig(i_lat,i_lon,:,1),[],1); a(isnan(a))=0;
    b=reshape(g.dens_sig(i_lat,i_lon,:,7),[],1); b(isnan(b))=0;
    fill([g.time fliplr(g.time)]', [a ; flipud(b)],[.9 .9 .9],'EdgeColor','none')
    a=reshape(g.dens_sig(i_lat,i_lon,:,2),[],1); a(isnan(a))=0;
    b=reshape(g.dens_sig(i_lat,i_lon,:,6),[],1); b(isnan(b))=0;
    fill([g.time fliplr(g.time)]', [a ; flipud(b)],[.8 .8 .8],'EdgeColor','none')
    a=reshape(g.dens_sig(i_lat,i_lon,:,3),[],1); a(isnan(a))=0;
    b=reshape(g.dens_sig(i_lat,i_lon,:,5),[],1); b(isnan(b))=0;
    fill([g.time fliplr(g.time)]', [a ; flipud(b)],[.7 .7 .7],'EdgeColor','none')
    plot(g.time,reshape(g.dens_sig(i_lat,i_lon,:,4),[],1),'k');
    xlabel('Date'); ylabel('Bird density'); axis tight; ylim([0 500]); box on; %set(gca,'YScale','log')

end



%% Figure 6: Cleaning
load('./data/d.mat');
load('./data/dc_corr.mat')

fig0=figure('position',[0 0 1000 400]); 

i_dc=find(strcmp({dc.name},'dedrs'));

for i_dc=1:numel(dc)
clf;
i_d=find(strcmp({d.name},dc(i_dc).name));
subplot(2,1,1); hold on;
imagesc(datenum(d(i_d).time), d(i_d).interval*(1/2:double(d(i_d).levels)), log(d(i_d).dens)','AlphaData',~isnan(d(i_d).dens'))
%plot([datenum(start_date) datenum(end_date)],[d(i_d).height d(i_d).height],'r','linewidth',2); caxis([-5 5])
xlabel('Date'); ylabel('Altitude [m]'); c=colorbar; c.Label.String='Bird Density [bird/m^3]'; caxis([0 5])
datetick('x'); axis([datenum(start_date) datenum(end_date-1) 0 5000]); set(gca, 'YDir', 'normal')

subplot(2,1,2); hold on;
imagesc(datenum(dc(i_dc).time), dc(i_dc).interval*(1/2:double(dc(i_dc).levels)), log(dc(i_dc).dens)', 'AlphaData',~isnan(dc(i_dc).dens)')
%plot([datenum(start_date) datenum(end_date)],[dc(i_d).height dc(i_d).height],'r','linewidth',2); caxis([-5 5])
xlabel('Date'); ylabel('Altitude [m]'); c=colorbar; c.Label.String='Bird Density [bird/m^3]'; caxis([0 5])
datetick('x'); axis([datenum(start_date) datenum(end_date-1) 0 5000]); set(gca, 'YDir', 'normal')
keyboard
end

% export_fig 'figure/paper/clenaning.eps' -eps













