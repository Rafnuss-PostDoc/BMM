%% Modeling Mird Migration as a Flow

% *Load data*

clear all;
load('data/Density_inference.mat','data');
load('data/Density_simulationMap_reassemble');
load('data/Density_estimationMap','g');
load('data/Flight_estimationMap','guv');
addpath('./functions/')
load coastlines;
col2=colormap(brewermap([],'Paired'));
%% 
% *Build Grid*

gext.lat = fillmissing([nan(1,1) ;g.lat ;nan(1,1)],'linear');
gext.lon = fillmissing([nan(1,1); g.lon; nan(1,1)],'linear');
gext.nlat = numel(gext.lat);
gext.nlon = numel(gext.lon);
[gext.lat2D, gext.lon2D] = ndgrid(gext.lat,gext.lon);

gext.time = sort( [g.time; g.day-.5]);

gext.nt = numel(gext.time);
gext.day_id=nan(gext.nt,1);
gext.day_id(ismember(gext.time,g.time)) = g.day_id;
gext.day_id(isnan(gext.day_id)) = gext.day_id(find(isnan(gext.day_id))+1);
tmp=nan(g.nlm,gext.nt);
tmp(:,ismember(gext.time,g.time)) = g.NNT;
gext.NNT=nan(gext.nlat,gext.nlon,gext.nt);
gext.NNT(padarray(repmat(g.latlonmask,1,1,gext.nt),[1 1 0], false))=tmp;
gext.mask_day = logical(gext.NNT>-1 & gext.NNT<1);
gext.time_b=nan(gext.nt,1);
gext.time_b(ismember(gext.time,g.time)) = g.time_b;
gext.latlonmask = padarray(g.latlonmask,[1 1],false);

dy = lldistkm([g.lat(1) g.lon(1)],[g.lat(2) g.lon(1)]);
dx = lldistkm([g.lat2D(:,1) g.lon2D(:,1)],[g.lat2D(:,1) g.lon2D(:,2)]);
area = repmat(dx*dy,1,g.nlon);
dxext = fillmissing([NaN;dx;NaN],'linear');
areaext = fillmissing(padarray(area,[1 1],nan),'linear',1);
areaext = fillmissing(areaext,'linear',2);

dt=15/60;

% *Compute flight speed and density*

vy = nan(g.nlat, g.nlon, gext.nt); vx=vy;
tmp = repmat(g.latlonmask,1,1,gext.nt); 
tmp(:,:,ismember(gext.time,g.day-.5))=false;
vy(tmp) = guv.v_est/1000*60*60; % m/s -> km/h (+) north, (-) south
vx(tmp) = guv.u_est/1000*60*60; % m/s -> km/h (+) east, (-) west

% * Transect definition*
transect_name={'North', 'East', 'Alpes', 'Spain' ,'Altlantic','England'};
transect(:,:,1)= [42 51 ; 43 82];
transect(:,:,2)= [20 45 ; 78 87];
transect(:,:,3)= [1 19 ; 40 84];
transect(:,:,4)= [1 9; 1 39];
transect(:,:,5) = [5 24; 1 18];
transect(:,:,6)= [25 43 ; 1 45];
sec_col={'r','b','y','g','c','k'};
cat=zeros(gext.nlat,gext.nlon);
for i_s=1:numel(transect_name)
    cat(transect(1,1,i_s):transect(1,2,i_s),transect(2,1,i_s):transect(2,2,i_s))=i_s;
end


%% * Compute for each realization*
nb_real = size(real_dens,4);

Ts_day_landing=nan(g.nat,nb_real);
Ts_day_takingoff=nan(g.nat,nb_real);
Ts_Fout_day_entering=nan(g.nat,nb_real);
Ts_Fout_day_leaving=nan(g.nat,nb_real);

Fout_sec_season=nan(numel(transect_name),4,nb_real);

for i_real = 1:nb_real
    tmp = real_dens(:,:,:,i_real);
    % remove for rain
    tmp(g.rain(repmat(g.latlonmask,1,1,g.nt))>data.mask_rain_thr) = 0;
    
    rho = nan(g.nlat, g.nlon, gext.nt);
    rho(:,:,~ismember(gext.time,g.day-.5)) = tmp;
    
    rho0 = rho;
    rho0(isnan(rho)) = 0; % bird/km^2

    % *Compute the Fluxes*

    Phiy_pad = padarray(rho .* vy .* repmat(dx,1,g.nlon,gext.nt) ,[1 0 0],nan); % bird/km^2 * km/h *km ->  bird/h
    Phix_pad = padarray(rho .* vx .* dy,[0 1 0],nan);

    Phiy_h = movmean(Phiy_pad,[0 1],1,'omitnan','Endpoints','discard'); 
    Phix_h = movmean(Phix_pad,[0 1],2,'omitnan','Endpoints','discard');

    Phiy_h_0=padarray(Phiy_h,[1 1 0],0);
    Phix_h_0=padarray(Phix_h,[1 1 0],0);

    Phiy_h_0(isnan(Phiy_h_0))=0;
    Phix_h_0(isnan(Phix_h_0))=0;

    dPhiy = diff(Phiy_h_0,1,1);
    dPhix = diff(Phix_h_0,1,2);

    F = (dPhiy + dPhix ).*dt; % bird/h * hr -> bird

    Fin = F; 
    Fin(~gext.mask_day)=0;
    Fin(repmat(~gext.latlonmask,1,1,gext.nt)) = 0;

    Fout = F;
    Fout(gext.mask_day)=0;
    Fout(repmat(gext.latlonmask,1,1,gext.nt)) = 0;

    W = diff(rho0 .* repmat(area,1,1,gext.nt),1,3) + Fin(2:end-1,2:end-1,1:end-1);

    
    % *Aggregate the Fluxes*
    
    % Separate the landing/takingoff
    W_landing=W; W_takingoff=W;
    W_landing(W_landing>=0)=0;
    W_takingoff(W_takingoff<=0)=0;
    
    % Separating entering and leaving
    Fout_leaving=Fout; Fout_entering=Fout;
    Fout_leaving(Fout_leaving>=0)=0;
    Fout_entering(Fout_entering<=0)=0;

    % Compute daily sum
    W_day_takingoff=nan(g.nlat,g.nlon,g.nat);
    W_day_landing=nan(g.nlat,g.nlon,g.nat);
    Fout_day_entering=nan(gext.nlat,gext.nlon,g.nat);
    Fout_day_leaving=nan(gext.nlat,gext.nlon,g.nat);
    for i_t=1:g.nat-1
        idt=gext.day_id(1:end-1)==i_t;
        W_day_takingoff(:,:,i_t)=nansum(W_takingoff(:,:,idt),3);
        W_day_landing(:,:,i_t)=nansum(W_landing(:,:,idt),3);
        Fout_day_entering(:,:,i_t)=nansum(Fout_entering(:,:,idt),3);
        Fout_day_leaving(:,:,i_t)=nansum(Fout_leaving(:,:,idt),3);
    end
    W_day_takingoff(W_day_takingoff==0)=nan;
    W_day_landing(W_day_landing==0)=nan;


    % *Compute timeseries*
    %Ts.W = reshape(nansum(nansum(W.W,1),2),1,[]);
    %Ts.landing = reshape(nansum(nansum(W_landing,1),2),1,[]);
    %Ts.takingoff = reshape(nansum(nansum(W_takingoff,1),2),1,[]);
    %Ts_Fout_entering = reshape(nansum(nansum(Fout_entering,1),2),1,[]);
    %Ts_Fout_leaving = reshape(nansum(nansum(Fout_leaving,1),2),1,[]);
    Ts_day_landing(:,i_real) = reshape(nansum(nansum(W_day_landing,1),2),1,[]);
    Ts_day_takingoff(:,i_real) = reshape(nansum(nansum(W_day_takingoff,1),2),1,[]);
    Ts_Fout_day_entering(:,i_real) = reshape(nansum(nansum(Fout_day_entering,1),2),1,[]);
    Ts_Fout_day_leaving(:,i_real) = reshape(nansum(nansum(Fout_day_leaving,1),2),1,[]);


    % * Fluxes per transect*
    for i_s=1:numel(transect_name)
        tmp = Fout;
        tmp(~(repmat(cat==i_s,1,1,gext.nt-1)))=0;
        Fout_sec(:,i_s) = reshape(nansum(nansum(tmp,1),2),1,[]);
    end
    
    % Flux per season
    tmp = Fout_sec(g.time<datetime('2018-07-01'),:); tmp(tmp>0)=0;
    Fout_sec_season(:,1,i_real)=nansum(tmp);
    tmp = Fout_sec(g.time<datetime('2018-07-01'),:); tmp(tmp<0)=0;
    Fout_sec_season(:,2,i_real)=nansum(tmp);
    tmp = Fout_sec(g.time>datetime('2018-07-01'),:); tmp(tmp>0)=0;
    Fout_sec_season(:,3,i_real)=nansum(tmp);
    tmp = Fout_sec(g.time>datetime('2018-07-01'),:); tmp(tmp<0)=0;
    Fout_sec_season(:,4,i_real)=nansum(tmp);

end

knownday = ismember(datenum(g.day), unique(data.day(data.day_id)));
Ts_day_landing(~knownday,:)=nan;
Ts_day_takingoff(~knownday,:)=nan;
Ts_Fout_day_leaving(~knownday,:)=nan;
Ts_Fout_day_entering(~knownday,:)=nan;  




%% 
% *Save*
%%
%save('data/SinkSource.mat','W','Fout','Ts','gext','vy','vx','rho');
%load('data/SinkSource.mat')

% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%% Results


%% Year-round accumulation of migratory birds on the ground
% Times serie daily

figure('position',[0 0 1000 1000]); subplot(3,1,[1 2]); box on; grid on
hold on; ylabel('Number of birdss (millions)')
fill(datenum([g.day(1); g.day; g.day(end)]),[0 cumsum((Ts.Fout.day.leaving+Ts.Fout.day.entering)/1000000,'omitnan') 0],'k','FaceAlpha',0.2);
h1=plot(datenum(g.day),cumsum((Ts_Fout_day_leaving+Ts_Fout_day_entering)/1000000,'omitnan'),'-k','linewidth',2 ,'DisplayName','Cumulative \Delta takingoff/landing = # on the ground');
h2=plot(datenum(g.day),Ts_day_landing/1000000,'Color',col2(6,:),'DisplayName','Landing');
h3=plot(datenum(g.day),Ts_day_takingoff/1000000,'Color',col2(2,:),'DisplayName','Takingoff');
h4=plot(datenum(g.day),Ts_Fout_day_leaving/1000000,'Color',col2(8,:),'DisplayName','Leaving');
h5=plot(datenum(g.day),Ts_Fout_day_entering/1000000,'Color',col2(10,:),'DisplayName','Entering +');
h6=plot(datenum(g.day),(Ts_Fout_day_leaving+Ts_Fout_day_entering)/1000000,'k','DisplayName','\Delta takingoff-landing = \Delta Leaving- Entering');
datetick('x'); axis tight; %legend([h2 h3 h4 h5 h6 h1]);

subplot(3,1,3); box on; grid on
hold on; ylabel('Number of birdss (millions)')
h1=fill(datenum([g.day ;flipud(g.day)]),-[cumsum((Ts.Fout.day.leaving)/1000000,'omitnan') fliplr(cumsum(-(Ts.Fout.day.entering)/1000000,'omitnan'))],'k','FaceAlpha',0.2,'DisplayName','Difference');
h2=plot(datenum(g.day),cumsum(-(Ts_Fout_day_leaving)/1000000,'omitnan'),'Color',col2(8,:),'linewidth',2 ,'DisplayName','Cumulative Leaving');
h3=plot(datenum(g.day),cumsum((Ts_Fout_day_entering)/1000000,'omitnan'),'Color',col2(10,:),'linewidth',2 ,'DisplayName','Cumulative Entering');
datetick('x'); xlabel('Date (2018)'); axis tight; %legend([h2 h3 h1]);

% Maximum Number of bird in a night
[max_takeoff,tmp] = max(Ts_day_takingoff(g.day<datetime('15-July-2018'),:));
median(max_takeoff/1000000)
median(g.day(tmp))
[max_takeoff,tmp] = max(Ts_day_takingoff);
median(max_takeoff/1000000)
median(g.day(tmp))

% How many night to reach 50% of all derpature
tmp = Ts_day_takingoff(g.day<datetime('15-July-2018') & ~isnan(Ts_day_takingoff'),:);
sum( (cumsum(sort(tmp,'descend'),'omitnan') ./ nansum(tmp)) < .5) 
tmp = Ts_day_takingoff(g.day>datetime('15-July-2018') & ~isnan(Ts_day_takingoff'),:);
sum( (cumsum(sort(tmp,'descend'),'omitnan') ./ nansum(tmp)) < .5) 
tmp = Ts_day_takingoff( ~isnan(Ts_day_takingoff'));
sum( (cumsum(sort(tmp,'descend'),'omitnan') ./ nansum(tmp)) < .5) 

% Accumulation of bird in Spring
spring_entering = nansum(Ts_Fout_day_entering(g.day<datetime('15-July-2018'),:)/1000000);
spring_leaving = nansum(Ts_Fout_day_leaving(g.day<datetime('15-July-2018'),:)/1000000);
spring_delta = spring_entering + spring_leaving;
median(spring_entering)
median(spring_leaving)
median(spring_delta)

% Accumulation (or dissipation) of bird in Autumns
autumn_entering = nansum(Ts_Fout_day_entering(g.day>datetime('15-July-2018'),:)/1000000);
autumn_leaving = nansum(Ts_Fout_day_leaving(g.day>datetime('15-July-2018'),:)/1000000);
autumn_delta = autumn_entering + autumn_leaving;
median(autumn_entering)
median(autumn_leaving)
median(autumn_delta)

% Recrutement = Flux out/Flux in = Bird in + reproduction / Bird in
-autumn_delta / spring_delta
%% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 

%% Seasonal Flow: In and Out



figure('position',[0 0 1000 400]); box on
bar(abs(Fout_sec_season)'/1000000,'stacked'); axis tight
xticklabels({'Spring -', 'Spring +','Autumn -','Autumn +'})
ylabel('Number of birds crossing the transect (-) outward, (+) inward [millions]')
legend(transect_name)
tmp = [sum(Fout_sec_season(:,1:2,:),2) sum(Fout_sec_season(:,3:4,:),2)]/1000000;

% Spring and Autum (col) fluxes over the transects
disp(round(mean(tmp,3)))

% Internal change
disp(round(sum(mean(tmp,3))))

% Sum of grouped transect
disp(round(mean([sum(tmp([1 2 end],:,:)) ; sum(tmp([3:5],:,:))],3)))
figure; subplot(1,2,1)
pie(flipud(abs(mean(tmp(:,1,:),3))),fliplr(transect_name))
subplot(1,2,2)
pie(flipud(abs(mean(tmp(:,2),3))),fliplr(transect_name))

%
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 