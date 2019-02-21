%% Map Estimation

%% 1. Load data
%
clear all; load('data/dc_corr.mat'); load coastlines; addpath('../functions/'); load('data/Flight_modelInf.mat'); 
load('data/Density_estimationMap.mat','g')

%% 6.3 Residu Kriging
%
Gneiting = @(dist,time,range_dist,range_time,delta,gamma,beta) 1./( (time./range_time).^(2.*delta) +1 ) .* exp(-(dist./range_dist).^(2.*gamma)./((time./range_time).^(2.*delta) +1).^(beta.*gamma) );

gtimenum=datenum(g.time);
uv_cov_parm_u=uv.cov.parm_u;
uv_cov_parm_v=uv.cov.parm_v;
g_mask_day=reshape(g.mask_day,g.nl,[]);
g_dateradar=g.dateradar;
ndc = numel(dc);

gDdist_sf = pdist2([[dc.lat]' [dc.lon]'], [g.lat2D(g.latlonmask) g.lon2D(g.latlonmask)], @lldistkm);
g_mask_day = g_mask_day(g.latlonmask(:),:)';

neigh_nb = 100;
dist_thr = 1000;%fminsearch(@(dist_thr) abs(.01-Gneiting(dist_thr,0, uv.cov.parm_u(3), uv.cov.parm_u(4), uv.cov.parm_u(5), uv.cov.parm_u(6), uv.cov.parm_u(7))) , uv.cov.parm_u(3) );
time_thr = 2;%fminsearch(@(time_thr) abs(.01-Gneiting(0,time_thr, uv.cov.parm_u(3), uv.cov.parm_u(4), uv.cov.parm_u(5), uv.cov.parm_u(6), uv.cov.parm_u(7))) , uv.cov.parm_u(4) );

guv_trans_est=cell(g.nat,1);
guv_trans_sig=guv_trans_est;
neighday=cell(g.nat,1);
sim=neighday;
uv_cov_C_u=cell(g.nat,1);
uv_cov_C_v=cell(g.nat,1);
uv_utrans=cell(g.nat,1);
uv_vtrans=cell(g.nat,1);
guv_utrans_est=cell(g.nat,1);
for i_d=1:g.nat
    neighday{i_d}  = find(ismember( data.dateradar, (i_d-1)*ndc+(1:ndc)));
    uv_cov_C_u{i_d} = uv.cov.C_u(neighday{i_d},neighday{i_d});
    uv_cov_C_v{i_d} = uv.cov.C_v(neighday{i_d},neighday{i_d});
    uv_utrans{i_d} =  uv.utrans(neighday{i_d});
    uv_vtrans{i_d} =  uv.vtrans(neighday{i_d});
    sim{i_d} = find(g_dateradar==i_d);
    guv_utrans_est{i_d} = nan(g.nlm,numel(sim{i_d}));
end
guv_utrans_sig = guv_utrans_est;
guv_vtrans_est = guv_utrans_est;
guv_vtrans_sig = guv_utrans_est;

% maxNumCompThreads('automatic');
% parpool(maxNumCompThreads());

parfor i_d=1:g.nat
   %maxNumCompThreads(1);
    for i_l=1:g.nlm
        tmpD = gDdist_sf(data.i_r(neighday{i_d}),i_l);
        neighloc = find(tmpD<dist_thr);
        tmpD=tmpD(neighloc);
        simloc = find(g_mask_day(sim{i_d},i_l));
        for i_sd=1:numel(simloc)
            tmpT = abs(data.time(neighday{i_d}(neighloc))-gtimenum(sim{i_d}(simloc(i_sd))));
            tmpidT = tmpT<time_thr;
            neighhour = neighloc(tmpidT);
            if isempty(neighhour)
                guv_utrans_est{i_d}(i_l,simloc(i_sd)) = 0;
                guv_utrans_sig{i_d}(i_l,simloc(i_sd)) = uv_cov_parm_u(1) + uv_cov_parm_u(2);
                guv_vtrans_est{i_d}(i_l,simloc(i_sd)) = 0;
                guv_vtrans_sig{i_d}(i_l,simloc(i_sd)) = uv_cov_parm_v(1) + uv_cov_parm_v(2);
            else
                Cab_u = uv_cov_parm_u(2).*Gneiting(tmpD(tmpidT), tmpT(tmpidT), uv_cov_parm_u(3), uv_cov_parm_u(4), uv_cov_parm_u(5), uv_cov_parm_u(6), uv_cov_parm_u(7));
                [~,id]=maxk(Cab_u,neigh_nb);
                lambda = uv_cov_C_u{i_d}(neighhour(id),neighhour(id))  \  Cab_u(id);
                guv_utrans_est{i_d}(i_l,simloc(i_sd)) = lambda' * uv_utrans{i_d}(neighhour(id));
                guv_utrans_sig{i_d}(i_l,simloc(i_sd)) = sqrt(uv_cov_parm_u(1) + uv_cov_parm_u(2) - lambda' * Cab_u(id));

                Cab_v = uv_cov_parm_v(2).*Gneiting(tmpD(tmpidT), tmpT(tmpidT), uv_cov_parm_v(3), uv_cov_parm_v(4), uv_cov_parm_v(5), uv_cov_parm_v(6), uv_cov_parm_v(7));
                [~,id]=maxk(Cab_v,neigh_nb);
                lambda = uv_cov_C_v{i_d}(neighhour(id),neighhour(id))  \  Cab_v(id);
                guv_vtrans_est{i_d}(i_l,simloc(i_sd)) = lambda' * uv_vtrans{i_d}(neighhour(id));
                guv_vtrans_sig{i_d}(i_l,simloc(i_sd)) = sqrt(uv_cov_parm_v(1) + uv_cov_parm_v(2) - lambda' * Cab_v(id));
            end
        end
    end
end

guv.utrans_est=nan(g.nlat,g.nlon,g.nt);
guv.vtrans_est=nan(g.nlat,g.nlon,g.nt);
guv.utrans_sig=nan(g.nlat,g.nlon,g.nt);
guv.vtrans_sig=nan(g.nlat,g.nlon,g.nt);

[latmask, lonmask]=ind2sub([g.nlat g.nlon],find(g.latlonmask));
for i_d=1:g.nat
    id = sub2ind([g.nlat g.nlon g.nt], repmat(latmask,numel(sim{i_d}),1), repmat(lonmask,numel(sim{i_d}),1), repelem(sim{i_d},g.nlm));
    guv.utrans_est(id) = guv_utrans_est{i_d}(:);
    guv.utrans_sig(id) = guv_utrans_sig{i_d}(:);
    guv.vtrans_est(id) = guv_vtrans_est{i_d}(:);
    guv.vtrans_sig(id) = guv_vtrans_sig{i_d}(:);
end

 
%% Back-transform
%
guv.u_est = guv.utrans_est *uv.trans.std(1) + uv.trans.mean(1);
guv.v_est = guv.vtrans_est *uv.trans.std(2) + uv.trans.mean(2);

guv.u_sig = sqrt(guv.utrans_sig.^2 * uv.trans.std(1)^2 );
guv.v_sig = sqrt(guv.vtrans_sig.^2 * uv.trans.std(2)^2 );

guv.u_q10 = guv.u_est .* norminv(.1).*guv.u_sig;
guv.v_q10 = guv.v_est .* norminv(.1).*guv.v_sig;
guv.u_q90 = guv.u_est .* norminv(.9).*guv.u_sig;
guv.v_q90 = guv.v_est .* norminv(.9).*guv.v_sig;


%%  Save
%
save('data/Flight_estimationMap','g','guv','-v7.3')
% load('data/Flight_estimationMap')