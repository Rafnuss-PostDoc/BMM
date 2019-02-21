%% Modeling Flight

%% 1. Import and create data
clear all; load('./data/dc_corr.mat'); load coastlines; addpath('../functions/'); warning('off')

%% 
% Convert the volumetric speed to a surface speed (|denss|)
for i_d=1:numel(dc)

    % Combine Speed and angle as complex number
    dc(i_d).u = dc(i_d).ff .* sind(dc(i_d).dd);
    dc(i_d).v = dc(i_d).ff .* cosd(dc(i_d).dd);

    % Aggregation according to elevation but weighted by the number of bird.
    % compute height bin size
    tmp=diff(double([0 max(0, (1:dc(i_d).levels).*dc(i_d).interval - dc(i_d).height )]))/1000;
    dc(i_d).denss = dc(i_d).dens*tmp';
    % Compute number of bird
    nb_bird=dc(i_d).dens .* repmat(tmp,size(dc(i_d).dens,1),1);
    % need to account only when speed and direction data are present
    nb_bird_2 = nb_bird;
    nb_bird_2(isnan(dc(i_d).ff) | isnan(dc(i_d).dd))=nan;
    % remove speed above 25.
    nb_bird_2(dc(i_d).ff>25)=nan;
    % check that at direction is representatife of at least 50% of the bird
    nb_bird_2(nansum(nb_bird_2,2)./nansum(nb_bird,2)<.5,:)=nan;
    % check that at direction is representatife of at least 3 altitude bins
    nb_bird_2(sum(nb_bird_2>0,2)<4,:)=nan;

    w_MTR = nb_bird_2 ./ repmat(nansum(nb_bird_2,2),1,size(nb_bird_2,2));
    
    dc(i_d).vs = nansum(w_MTR.* dc(i_d).v,2);
    dc(i_d).us = nansum(w_MTR.* dc(i_d).u,2);
    dc(i_d).scoret = datenum(dc(i_d).time'-mean([dc(i_d).sunrise;dc(i_d).sunset])) ./ datenum(dc(i_d).sunrise-dc(i_d).sunset)*2;
end

%% 
% Create the data structure |data|
dens_lim = 0;
data=[];
uniqueDate = round(datenum(datetime(start_date-1:end_date-1)));
for i_d=1:numel(dc)
    id = dc(i_d).scoret'>=-1 & dc(i_d).scoret'<=1 & dc(i_d).denss>dens_lim & dc(i_d).vs~=0;
    d=table();
    d.v=dc(i_d).vs(id);
    d.u=dc(i_d).us(id);
    d.dens=dc(i_d).denss(id);
    d.scoret = dc(i_d).scoret(id)';
    d.i_r=i_d*ones(sum(id),1);
    d.lat=dc(i_d).lat*ones(sum(id),1);
    d.lon=dc(i_d).lon*ones(sum(id),1);
    d.time = datenum(dc(i_d).time(id));
    tmp = round(datenum(dc(i_d).sunrise(id)))';
    [~,id]=ismember(tmp,uniqueDate);
    d.dateradar = sub2ind([numel(dc) numel(uniqueDate)],repmat(i_d,numel(tmp),1),id);
    data = [data ; d ];
end

clear d dens_lim 


%% Normalization
%

figure;
subplot(2,1,1);histogram(data.v)
subplot(2,1,2);histogram(data.u)

uv.trans.mean = [mean(data.v) mean(data.u)];
uv.trans.std = [std(data.v) std(data.u)];

uv.vtrans = ( data.v-uv.trans.mean(1)) / uv.trans.std(1);
uv.utrans = ( data.u-uv.trans.mean(2)) / uv.trans.std(2);

figure; 
subplot(2,1,1); hold on;
histogram(uv.vtrans,'Normalization','pdf')
x=-3:0.1:3;
plot(x,normpdf(x))
subplot(2,1,2); hold on
histogram(uv.utrans,'Normalization','pdf')
plot(x,normpdf(x))



%% 5.2 Built matrix of distance
% Built the matrix of distances of the data. Distence (or difference) in lattitude, 
% longitude, time and value. 

Ddist_sf=squareform(pdist(repmat([[dc.lat]' [dc.lon]'], numel(uniqueDate),1),@lldistkm));

uv.D=table();
for i_d=1:numel(uniqueDate)
    d=table();
    id  = find(ismember( data.dateradar, sub2ind([numel(dc), numel(uniqueDate)],1:numel(dc),i_d*ones(1,numel(dc)))));
    d.Ddist = reshape(Ddist_sf(data.i_r(id),data.i_r(id)'),[],1);
    d.Dtime = reshape(squareform(pdist(data.time(id))),[],1);
    d.Did1 = repmat(id,numel(id),[]);
    d.Did2 = repelem(id,numel(id));
    
    A = ones(numel(id));
    B = tril(A);
    uv.D=[uv.D ; d(B(:)==1,:)];
end


%% 5.3 Empirical Variogram
%
uv.cov.d=[0 eps 200 300 400 600 1000 1500];
uv.cov.t=[0 eps .015 .025 .045 .08 .1 .2 .25  .35]; 
[uv.cov.D,uv.cov.T]=meshgrid(uv.cov.d(1:end-1)+diff(uv.cov.d)/2,uv.cov.t(1:end-1)+diff(uv.cov.t)/2);

uv.cov.emp_grid_u=nan(size(uv.cov.D));
uv.cov.emp_grid_v=nan(size(uv.cov.D));
nb=nan(size(uv.cov.D));
for i_d=1:numel(uv.cov.d)-1
    id1 = find(uv.D.Ddist>=uv.cov.d(i_d) & uv.D.Ddist<uv.cov.d(i_d+1));
    for i_t=1:numel(uv.cov.t)-1
        id2 = id1(uv.D.Dtime(id1)>=uv.cov.t(i_t) & uv.D.Dtime(id1)<uv.cov.t(i_t+1));
        uv.cov.emp_grid_u(i_t,i_d)=mean(uv.utrans(uv.D.Did1(id2)) .* uv.utrans(uv.D.Did2(id2)));
        uv.cov.emp_grid_v(i_t,i_d)=mean(uv.vtrans(uv.D.Did1(id2)) .* uv.vtrans(uv.D.Did2(id2)));
        nb(i_t,i_d)=numel(id2);
    end
end

figure; 
subplot(1,2,1); surf(uv.cov.D,uv.cov.T,uv.cov.emp_grid_u);  xlabel('Distance [km]'); ylabel('Time [Days]');
subplot(1,2,2); surf(uv.cov.D,uv.cov.T,uv.cov.emp_grid_v);  xlabel('Distance [km]'); ylabel('Time [Days]');
%subplot(1,3,3); surf(uv.cov.D,uv.cov.T,nb);  xlabel('Distance [km]'); ylabel('Time [Days]');


%% 5.4 Calibrate ampli.cov.parmriance 
%
Gneiting = @(dist,time,range_dist,range_time,delta,gamma,beta) 1./( (time./range_time).^(2.*delta) +1 ) .* exp(-(dist./range_dist).^(2.*gamma) );
Gneiting_fit = @(parm) parm(2).*Gneiting(uv.cov.D(:),uv.cov.T(:),parm(3),parm(4),parm(5),parm(6),parm(7));

parm0 =[.1 .9 300 5 .5 .5 .5];
parm_min=[0 0 0.0001 0 0 0 0 ]';
parm_max=[nancov(uv.utrans) nancov(uv.utrans) 3000 10 1 1 1]';

rmse_u = @(parm) sum( ( Gneiting_fit(parm) -  uv.cov.emp_grid_u(:) ).^2 );
uv.cov.parm_u = fmincon(rmse_u,parm0,[],[],[],[],parm_min,parm_max,[]); % ,optimset('PlotFcn',@optimplotfval)
uv.cov.parm_u(1) = nancov(uv.utrans)-uv.cov.parm_u(2);

parm_max=[nancov(uv.vtrans) nancov(uv.vtrans) 3000 10 1 1 1]';
rmse_u = @(parm) sum( ( Gneiting_fit(parm) -  uv.cov.emp_grid_v(:) ).^2 );
uv.cov.parm_v = fmincon(rmse_u,parm0,[],[],[],[],parm_min,parm_max,[]); % ,optimset('PlotFcn',@optimplotfval)
uv.cov.parm_v(1) = nancov(uv.vtrans)-uv.cov.parm_v(2);

sprintf('%.3f ',uv.cov.parm_u)
sprintf('%.3f ',uv.cov.parm_v)


% Figure
[tmpD,tmpT] = meshgrid(uv.cov.d(1):1:uv.cov.d(end),uv.cov.t(1):0.01:uv.cov.t(end));
figure; 

subplot(1,2,1); hold on;
tmpCOV=uv.cov.parm_u(2).*Gneiting(tmpD(:), tmpT(:), uv.cov.parm_u(3), uv.cov.parm_u(4), uv.cov.parm_u(5), uv.cov.parm_u(6), uv.cov.parm_u(7));
tmpCOV(1)=tmpCOV(1)+uv.cov.parm_u(1);
tmpCOV=reshape(tmpCOV,size(tmpD,1),size(tmpD,2));
s=surf(tmpD,tmpT,tmpCOV);
surf(uv.cov.D,uv.cov.T,uv.cov.emp_grid_u);  xlabel('Distance [km]'); ylabel('Time [Days]');
s.EdgeColor='none'; view(3)

subplot(1,2,2); hold on;
tmpCOV=uv.cov.parm_v(2).*Gneiting(tmpD(:), tmpT(:), uv.cov.parm_v(3), uv.cov.parm_v(4), uv.cov.parm_v(5), uv.cov.parm_v(6), uv.cov.parm_v(7));
tmpCOV(1)=tmpCOV(1)+uv.cov.parm_v(1);
tmpCOV=reshape(tmpCOV,size(tmpD,1),size(tmpD,2));
s=surf(tmpD,tmpT,tmpCOV);
surf(uv.cov.D,uv.cov.T,uv.cov.emp_grid_v);  xlabel('Distance [km]'); ylabel('Time [Days]');
s.EdgeColor='none'; view(3)


%% 5.5 Kriging residual
%%
tmp = uv.cov.parm_u(2).*Gneiting(uv.D.Ddist, uv.D.Dtime, uv.cov.parm_u(3), uv.cov.parm_u(4), uv.cov.parm_u(5), uv.cov.parm_u(6), uv.cov.parm_u(7));
tmp = sparse(uv.D.Did1,uv.D.Did2,tmp);
uv.cov.C_u=(tmp'+tmp);
uv.cov.C_u(1:size(tmp,1)+1:end) = uv.cov.parm_u(1)+uv.cov.parm_u(2);

tmp = uv.cov.parm_v(2).*Gneiting(uv.D.Ddist, uv.D.Dtime, uv.cov.parm_v(3), uv.cov.parm_v(4), uv.cov.parm_v(5), uv.cov.parm_v(6), uv.cov.parm_v(7));
tmp = sparse(uv.D.Did1,uv.D.Did2,tmp);
uv.cov.C_v=(tmp'+tmp);
uv.cov.C_v(1:size(tmp,1)+1:end) = uv.cov.parm_v(1)+uv.cov.parm_v(2);

uv.utrans_est=nan(size(uv.utrans));
uv.utrans_sig=nan(size(uv.utrans));
uv.vtrans_est=nan(size(uv.vtrans));
uv.vtrans_sig=nan(size(uv.vtrans));

k=100;

for i_d=1:numel(uniqueDate)
    neigh  = find(ismember( data.dateradar, sub2ind([numel(dc), numel(uniqueDate)],1:numel(dc),i_d*ones(1,numel(dc)))));

    for i_s=1:numel(neigh)
        
        neigh_s = neigh(data.dateradar(neigh)~=data.dateradar(neigh(i_s)));
        Cab = full(uv.cov.C_u(neigh_s,neigh(i_s)));
        [~,id]=maxk(Cab,k);
        neigh_s=neigh_s(id);

        % Ordinary
        lambdaU = [uv.cov.C_u(neigh_s,neigh_s) ones(numel(neigh_s),1); ones(1,numel(neigh_s)) 0] \ [ uv.cov.C_u(neigh_s,neigh(i_s)) ; 1];
        uv.utrans_est(neigh(i_s)) = lambdaU(1:end-1)' * uv.utrans(neigh_s);
        uv.utrans_sig(neigh(i_s)) = sqrt( uv.cov.parm_u(1) + uv.cov.parm_u(2) - lambdaU' * [uv.cov.C_u(neigh_s,neigh(i_s)); 1] );
       
        lambdaV = [uv.cov.C_v(neigh_s,neigh_s) ones(numel(neigh_s),1); ones(1,numel(neigh_s)) 0] \ [ uv.cov.C_v(neigh_s,neigh(i_s)) ; 1];
        uv.vtrans_est(neigh(i_s)) = lambdaV(1:end-1)' * uv.vtrans(neigh_s);
        uv.vtrans_sig(neigh(i_s)) = sqrt( uv.cov.parm_v(1) + uv.cov.parm_v(2) - lambdaV' * [uv.cov.C_v(neigh_s,neigh(i_s)); 1] );
    end
end


%% 
% Figure

figure('position',[0 0 1000 400]); 
subplot(2,1,1); hold on;
plot(uv.utrans,'k');
plot(uv.utrans_est,'r');
plot(uv.utrans_est-uv.utrans_sig,'--r')
plot(uv.utrans_est+uv.utrans_sig,'--r')

subplot(2,1,2); hold on;
plot(uv.vtrans,'k');
plot(uv.vtrans_est,'r');
plot(uv.vtrans_est-uv.vtrans_sig,'--r')
plot(uv.vtrans_est+uv.vtrans_sig,'--r')


err_norm_u = (uv.utrans_est-uv.utrans)./(uv.utrans_sig);
err_norm_v = (uv.vtrans_est-uv.vtrans)./(uv.vtrans_sig);

figure('position',[0 0 1000 400]); hold on;  histogram( err_norm_u ); histogram( err_norm_v );
legend(['Normalized error of kriging: mean=' num2str(nanmean(err_norm_u)) ' and std=' num2str(nanstd(err_norm_u))], ['Normalized error of kriging: mean=' num2str(nanmean(err_norm_v)) ' and std=' num2str(nanstd(err_norm_v))]);


% Map of normalized error of estimation per radar
figure('position',[0 0 1000 400]); hold on
worldmap([min([dc.lat]) max([dc.lat])], [min([dc.lon]) max([dc.lon])]); 
plotm(coastlat, coastlon,'k');
geoshow('worldrivers.shp','Color', 'blue')
scatterm([dc.lat],[dc.lon],splitapply(@std,err_norm_u,data.i_r)*200,splitapply(@mean,err_norm_u,data.i_r),'filled','MarkerEdgeColor','k'); 
scatterm(min([dc.lat])*[1 1 1 1],min([dc.lon])*[1 1 1 1],[1.5 1 0.5 0.25]*200,'filled','MarkerEdgeColor','k'); 
colorbar;

% figure('position',[0 0 1000 400]); histogram( err_norm );
% legend(['Normalized error of kriging: mean=' num2str(mean(err_norm)) ' and std=' num2str(std(err_norm))]);
%     
% figure('position',[0 0 1000 400]); histogram( err_norm );
% boxplot(err_norm,data.i_r) 

%% 5.6 Back Normalized
%

uv.u_est = uv.utrans_est *uv.trans.std(1) + uv.trans.mean(1);
uv.v_est = uv.vtrans_est *uv.trans.std(2) + uv.trans.mean(2);

uv.u_q10 = uv.utrans_est *norminv(.1)*uv.trans.std(1) + uv.trans.mean(1);
uv.v_q10 = uv.vtrans_est *norminv(.1)*uv.trans.std(2) + uv.trans.mean(2);

uv.u_q90 = uv.utrans_est *norminv(.9)*uv.trans.std(1) + uv.trans.mean(1);
uv.v_q90 = uv.vtrans_est *norminv(.9)*uv.trans.std(2) + uv.trans.mean(2);

figure('position',[0 0 1000 400]); 
subplot(2,1,1); hold on;
plot(data.u,'k');
plot(uv.u_est,'r');
plot(uv.u_q10,'--r')
plot(uv.u_q90,'--r')

subplot(2,1,2); hold on;
plot(data.v,'k');
plot(uv.v_est,'r');
plot(uv.v_q10,'--r')
plot(uv.v_q90,'--r')


%% 7. Save
%
save('data/Flight_modelInf.mat','data','uv','-v7.3')
% load('data/Flight_modelInf.mat')