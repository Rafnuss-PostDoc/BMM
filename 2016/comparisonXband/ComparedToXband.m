%% Comparaison to Xband
% 1. Setup

clear all; addpath('functions/'); load('data/dc_corr.mat');

% Set-up grid

pt.lat=[50.884889, 47.126356];
pt.lon=[2.544694, 8.193466];
pt.time = start_date:1/24/4:end_date;
pt.atime=unique(round(datenum(pt.time)));

pt.nl=numel(pt.lat);
pt.nt=numel(pt.time);
pt.nat=numel(pt.atime);

pt.time2D=repmat(pt.time,2,1);

% Sunrise and sunset

load('./data/sunrisesunset_grid.mat')

pt.sunrise = nan(pt.nl,pt.nt);
pt.sunset = pt.sunrise;
for i_t=2:numel(tim_d)
    id = pt.time>tim_d(i_t)-1 & pt.time<=tim_d(i_t);
    
    F=griddedInterpolant({lat_d,lon_d},sunrs_e(:,:,i_t-1));
    pt.sunset(:,id)=repmat(F(pt.lat,pt.lon),1,1,sum(id));
    
    F=griddedInterpolant({lat_d,lon_d},sunrs_b(:,:,i_t));
    pt.sunrise(:,id)=repmat(F(pt.lat,pt.lon),1,1,sum(id));
end

% Score time of night

pt.scoret = (datenum(pt.time2D)-mean(cat(3,pt.sunrise,pt.sunset),3)) ./ (pt.sunrise-pt.sunset)*2;
pt.scoret(pt.scoret<-1.1)=-1.1;
pt.scoret(pt.scoret>1.1)=1.1;

% Mask of the day

pt.mask_day =  pt.scoret>-1 & pt.scoret<1;

% Group g.time by night

[~,pt.dateradar] = ismember(round(datenum(pt.sunrise(1,:))),pt.atime);
pt.dateradar(pt.dateradar==0)=1;

% Clear the masking variable

clear sunrs_b sunrs_e tim_d lat_d lon_d F 


%% 2. Density Estimation
% load data

load('data/Density_modelInf.mat'); addpath('../functions/')

% Amplitude Kriging
% Built ampli.cov.parmriance of the radar data
ptampli.Ddist_sf = pdist2(repmat([[dc.lat]' [dc.lon]'], pt.nat,1), [pt.lat(:) pt.lon(:)], @lldistkm);
ptampli.Dtime_sf = pdist2(repelem(pt.atime',numel(dc),1), pt.atime');
ptampli.Dtime_sf = ptampli.Dtime_sf(~ampli.isnanA,:);
ptampli.Ddist_sf = ptampli.Ddist_sf(~ampli.isnanA,:);

% Initialized the simulated grid
ptampli.An_est=nan(pt.nat,pt.nl);
ptampli.An_sig=ptampli.An_est;

% Perform the kriging

Gneiting = @(dist,time,range_dist,range_time,delta,gamma,beta) 1./( (time./range_time).^(2.*delta) +1 ) .* exp(-(dist./range_dist).^(2.*gamma)./((time./range_time).^(2.*delta) +1).^(beta.*gamma) );
wradius=4;

ampli_AnNan = ampli.An(~isnan(ampli.An));

for i_l=1:pt.nl
    for i_t=1:pt.nat
        id2=ptampli.Ddist_sf(:,i_l)<ampli.cov.parm(3)*wradius & ptampli.Dtime_sf(:,i_t)<ampli.cov.parm(4)*wradius;
        Cab = ampli.cov.parm(2).*Gneiting(ptampli.Ddist_sf(id2,i_l), ptampli.Dtime_sf(id2,i_t), ampli.cov.parm(3), ampli.cov.parm(4), ampli.cov.parm(5), ampli.cov.parm(6), ampli.cov.parm(7));
        lambda = ampli.cov.C(id2,id2)  \  Cab;
        ptampli.An_est(i_t,i_l) = lambda' * ampli_AnNan(id2);
        ptampli.An_sig(i_t,i_l) = ampli.cov.parm(1) + ampli.cov.parm(2) - lambda' * Cab;
    end
end

ptampli.A_est = ptampli.An_est*ampli.ns.std;
ptampli.A_sig = sqrt(ptampli.An_sig.^2 * ampli.ns.std.^2);

% Residu Kriging

Gneiting = @(dist,time,range_dist,range_time,delta,gamma,beta) 1./( (time./range_time).^(2.*delta) +1 ) .* exp(-(dist./range_dist).^(2.*gamma)./((time./range_time).^(2.*delta) +1).^(beta.*gamma) );
pttimenum=datenum(pt.time);
ptDdist_sf = pdist2([[dc.lat]' [dc.lon]'], [pt.lat(:) pt.lon(:)], @lldistkm);

neighday=cell(pt.nat,1);
sim=cell(pt.nat,1);
res_cov_C=cell(pt.nat,1);
res_rn=cell(pt.nat,1);
for i_d=1:pt.nat
    neighday{i_d}  = find(ismember( data.dateradar, (i_d-1)* numel(dc)+(1:numel(dc))));
    res_cov_C{i_d} = res.cov.C(neighday{i_d},neighday{i_d});
    res_rn{i_d} =  res.rn(neighday{i_d});
    sim{i_d} = find(pt.dateradar==i_d);
end

ptres.rn_est=nan(pt.nl,pt.nt);
ptres.rn_sig=ptres.rn_est;

for i_d=1:pt.nat
    for i_l=1:pt.nl
        tmpD = ptDdist_sf(data.i_r(neighday{i_d}),i_l);
        neighloc = find(tmpD<res.cov.parm(3)*wradius);
        tmpD=tmpD(neighloc);
        simloc = find(pt.mask_day(i_l,sim{i_d}));
        for i_sd=1:numel(simloc)
            tmpT = abs(data.time(neighday{i_d}(neighloc))-pttimenum(sim{i_d}(simloc(i_sd))));
            tmpidT = tmpT<res.cov.parm(4)*wradius;
            neighhour = neighloc(tmpidT);
            if isempty(neighhour)
                ptres.rn_est(i_l,sim{i_d}(simloc(i_sd))) = 0;
                ptres.rn_sig(i_l,sim{i_d}(simloc(i_sd))) = res.cov.parm(1) + res.cov.parm(2);
            else
                Cab = res.cov.parm(2).*Gneiting(tmpD(tmpidT), tmpT(tmpidT), res.cov.parm(3), res.cov.parm(4), res.cov.parm(5), res.cov.parm(6), res.cov.parm(7));
                lambda = res_cov_C{i_d}(neighhour,neighhour)  \  Cab;
                ptres.rn_est(i_l,sim{i_d}(simloc(i_sd))) = lambda' * res_rn{i_d}(neighhour);
                ptres.rn_sig(i_l,sim{i_d}(simloc(i_sd))) = res.cov.parm(1) + res.cov.parm(2) - lambda' * Cab;
            end
        end
    end
end
 
% Back-transform
ptres.r_est = ptres.rn_est.*res.ns.std(pt.scoret);
ptres.r_sig = sqrt(ptres.rn_sig.^2.*res.ns.std(pt.scoret).^2);


% Reassemble
t = pt.lat.*trend.p(1)+trend.p(2);
t = repmat(t',1,pt.nt);
A = ptampli.A_est(pt.dateradar,:)';
p = polyval(curve.p,pt.scoret);

pt.denstrans_est = t + A + p + ptres.r_est;
pt.dens_est = (pt.denstrans_est).^(1/pow_a);

pt.denstrans_sig = sqrt(ptampli.A_sig(pt.dateradar,:)'.^2 + ptres.r_sig.^2);


pt.dens_q10 = ( pt.denstrans_est+norminv(.1).*pt.denstrans_sig ).^(1/pow_a);
pt.dens_q90 = ( pt.denstrans_est+norminv(.9).*pt.denstrans_sig ).^(1/pow_a);


nsig=7;
pt.dens_sig = (repmat(pt.denstrans_est,1,1,nsig) + reshape(pt.denstrans_sig(:)*(-3:3),pt.nl,pt.nt,nsig)).^(1/pow_a);


% Figure
figure('position',[0 0 1000 600]); hold on;
subplot(2,1,1); xlabel('Time (day)'); ylabel('Bird/km^2')
plot(pt.time, reshape(pt.dens_sig(1,:,:),pt.nt,[]));
xlim([start_date end_date])
title('Herzeele')
subplot(2,1,2); xlabel('Time (day)'); ylabel('Bird/km^2')
plot(pt.time, reshape(pt.dens_sig(2,:,:),pt.nt,[]));
xlim([start_date end_date])
title('Sempach')

% Save
save('data/Density_validationRadars','pt','ptampli','ptres','-v7.3')


%% 1. Flight Speed Estimation at the X-band radar location

% Import Estimation data
load('data/Density_validationRadars.mat');
load('data/Density_modelInf.mat')
 
% Import x-band data
T_Herzeele = readtable(['comparisonXband/MTR_SorL_allBirds_hourly_Hrange50_1500asl_Hstep1450_inclFlight_13dBadj.csv'],'TreatAsEmpty','NA');
T_Herzeele(strcmp(T_Herzeele.OperMod,'L'),:)=[];
T_Herzeele.dens = T_Herzeele.mtr/60/60./(T_Herzeele.meanSpeed./1000);
T_Herzeele.denstrans = T_Herzeele.dens.^pow_a;

T_Sem = readtable(['comparisonXband/SEM_MTR_allBirds_SeptOct2016_1h_2000m_inclFlight_sepPulse.csv'],'TreatAsEmpty','NA');
T_Sem(strcmp(T_Sem.oper_mod,'L'),:)=[];
T_Sem(strcmp(T_Sem.oper_mod,'M'),:)=[];
T_Sem.dens = T_Sem.mtr/60/60./(T_Sem.meanSpeed./1000);
T_Sem.denstrans = T_Sem.dens.^pow_a;



figure('position',[0 0 1000 600]); 
subplot(2,1,1);hold on;
plot(pt.time, pt.denstrans_est(1,:) -3*pt.denstrans_sig(1,:),'--','Color',0.8*[1 1 1]);
plot(pt.time, pt.denstrans_est(1,:) - 2*pt.denstrans_sig(1,:),'--','Color',0.4*[1 1 1]);
plot(pt.time, pt.denstrans_est(1,:) - pt.denstrans_sig(1,:) ,'--','Color',0.2*[1 1 1]);
plot(pt.time, pt.denstrans_est(1,:),'-', 'Color',0*[1 1 1]);
plot(pt.time, pt.denstrans_est(1,:) + pt.denstrans_sig(1,:),'--','Color',0.2*[1 1 1]);
plot(pt.time, pt.denstrans_est(1,:) + 2*pt.denstrans_sig(1,:),'--','Color',0.4*[1 1 1]);
plot(pt.time, pt.denstrans_est(1,:) + 3*pt.denstrans_sig(1,:),'--','Color',0.8*[1 1 1]);
plot(T_Herzeele.Tint,T_Herzeele.denstrans,'-r','linewidth',2);
xlim([start_date end_date]);xlabel('Time (day)'); ylabel('Z^p')
subplot(2,1,2);hold on;
plot(pt.time, pt.denstrans_est(2,:) -3*pt.denstrans_sig(2,:),'--','Color',0.8*[1 1 1]);
plot(pt.time, pt.denstrans_est(2,:) - 2*pt.denstrans_sig(2,:),'--','Color',0.4*[1 1 1]);
plot(pt.time, pt.denstrans_est(2,:) - pt.denstrans_sig(2,:) ,'--','Color',0.2*[1 1 1]);
plot(pt.time, pt.denstrans_est(2,:),'-', 'Color',0*[1 1 1]);
plot(pt.time, pt.denstrans_est(2,:) + pt.denstrans_sig(2,:),'--','Color',0.2*[1 1 1]);
plot(pt.time, pt.denstrans_est(2,:) + 2*pt.denstrans_sig(2,:),'--','Color',0.4*[1 1 1]);
plot(pt.time, pt.denstrans_est(2,:) + 3*pt.denstrans_sig(2,:),'--','Color',0.8*[1 1 1]);
plot(T_Sem.Tint,T_Sem.denstrans,'-r','linewidth',2);
xlim([start_date end_date]);xlabel('Time (day)'); ylabel('Z^p')


figure('position',[0 0 1000 600]); hold on;
subplot(2,1,1);hold on;
plot(pt.time, reshape(pt.dens_sig(1,:,1),pt.nt,[]),'--','Color',0.8*[1 1 1]);
plot(pt.time, reshape(pt.dens_sig(1,:,2),pt.nt,[]),'--','Color',0.4*[1 1 1]);
plot(pt.time, reshape(pt.dens_sig(1,:,3),pt.nt,[]),'--','Color',0.2*[1 1 1]);
plot(pt.time, reshape(pt.dens_sig(1,:,4),pt.nt,[]),'-', 'Color',0*[1 1 1]);
plot(pt.time, reshape(pt.dens_sig(1,:,5),pt.nt,[]),'--','Color',0.2*[1 1 1]);
plot(pt.time, reshape(pt.dens_sig(1,:,6),pt.nt,[]),'--','Color',0.4*[1 1 1]);
plot(pt.time, reshape(pt.dens_sig(1,:,7),pt.nt,[]),'--','Color',0.8*[1 1 1]);
plot(T_Herzeele.Tint,T_Herzeele.dens,'-r','linewidth',2);
xlim([start_date end_date]);xlabel('Time (day)'); ylabel('Bird/km^2')
subplot(2,1,2);hold on;
plot(pt.time, reshape(pt.dens_sig(2,:,1),pt.nt,[]),'--','Color',0.8*[1 1 1]);
plot(pt.time, reshape(pt.dens_sig(2,:,2),pt.nt,[]),'--','Color',0.4*[1 1 1]);
plot(pt.time, reshape(pt.dens_sig(2,:,3),pt.nt,[]),'--','Color',0.2*[1 1 1]);
plot(pt.time, reshape(pt.dens_sig(2,:,4),pt.nt,[]),'-', 'Color',0*[1 1 1]);
plot(pt.time, reshape(pt.dens_sig(2,:,5),pt.nt,[]),'--','Color',0.2*[1 1 1]);
plot(pt.time, reshape(pt.dens_sig(2,:,6),pt.nt,[]),'--','Color',0.4*[1 1 1]);
plot(pt.time, reshape(pt.dens_sig(2,:,7),pt.nt,[]),'--','Color',0.8*[1 1 1]);
plot(T_Sem.Tint,T_Sem.dens,'-r','linewidth',2);
xlim([start_date end_date]);xlabel('Time (day)'); ylabel('Bird/km^2')
%i_d=strcmp({dc.name},'bejab');
%plot(datetime(data.time(data.i_r==i_d),'convertfrom','datenum'),data.dens(data.i_r==i_d),'g')

% Score

[~,idx] = ismember(T_Herzeele.Tint,pt.time);
idx2=idx>0;
err_norm_her = (pt.denstrans_est(1,idx(idx2))' - T_Herzeele.denstrans(idx2)) ./ pt.denstrans_sig(1,idx(idx2))';

[~,idx] = ismember(T_Sem.Tint,pt.time);
idx2=idx>0;
err_norm_sem = (pt.denstrans_est(2,idx(idx2))' - T_Sem.denstrans(idx2)) ./ pt.denstrans_sig(2,idx(idx2))';

figure('position',[0 0 1000 400]); hold on; histogram( err_norm_her ); histogram( err_norm_sem );
legend(['Normalized error of kriging: mean=' num2str(nanmean(err_norm_her)) ' and std=' num2str(nanvar(err_norm_her))],['Normalized error of kriging: mean=' num2str(nanmean(err_norm_sem)) ' and std=' num2str(nanvar(err_norm_sem))]);



