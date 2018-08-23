    
%% Clear and load
clear all; load coastlines; addpath('functions/'); 
load('data/dc_corr'); 
load('data/krigingAmplitudeResidus')
load('data/krigingAmplitudeResidusGrid')

%% Amplitude
wradius=4;
neighg_nb=100;
neighr_nb=50;
Gneiting = @(dist,time,range_dist,range_time,delta,gamma,beta) 1./( (time./range_time).^(2.*delta) +1 ) .* exp(-(dist./range_dist).^(2.*gamma)./((time./range_time).^(2.*delta) +1).^(beta.*gamma) );

% 1. Creation of the grid and path
[LAT,LON,AT]= ndgrid(1:g.nlat,1:g.nlon,1:g.nat);
tmp=Inf*ones(g.nlat,g.nlon);
tmp(g.mask_water&g.mask_distrad)=nan;
Path=repmat(tmp,1,1,g.nat);

slat = 1:ceil(log(g.nlat+1)/log(2));
slon = 1:ceil(log(g.nlon+1)/log(2));
sat = 1:ceil(log(g.nat+1)/log(2));
sn = max([numel(slat), numel(slon) numel(sat)]);
nb = nan(sn,1);
start = zeros(sn+1,1);
path = nan(sum(~~(Path(:)==0)),1);
ds = 2.^(sn-1:-1:0);
for i_scale = 1:sn
   [LAT_s,LON_s,AT_s] = ndgrid(1:ds(i_scale):g.nlat,1:ds(i_scale):g.nlon,1:ds(i_scale):g.nat); % matrix coordinate
   id = find(isnan(Path(:)) & ismember([LAT(:) LON(:) AT(:)], [LAT_s(:) LON_s(:) AT_s(:)], 'rows'));
   nb(i_scale) = numel(id);
   start(i_scale+1) = start(i_scale)+nb(i_scale);
   path( start(i_scale)+(1:nb(i_scale)) ) = id(randperm(nb(i_scale)));
   Path(path( start(i_scale)+(1:nb(i_scale)) )) = start(i_scale)+(1:nb(i_scale));
end

% Ditance Matrix
Ddist_gg = squareform(pdist([g.lat2D(:) g.lon2D(:)], @lldistkm));
Ddist_gr = pdist2(repmat([[dc.lat]' [dc.lon]'], g.nat,1), [g.lat2D(:) g.lon2D(:)], @lldistkm);
Ddist_gr = Ddist_gr(~ampli.isnanA,:);
Dtime_gg = squareform(pdist(g.atime'));
Dtime_gr = pdist2(repelem(g.atime',numel(dc),1), g.atime');
Dtime_gr = Dtime_gr(~ampli.isnanA,:);

nan(22176,22176,7);
ntmax=ceil(ampli.cov.parm(4)*wradius);
Cgp_i=ampli.cov.parm(2).*Gneiting(repmat(Ddist_gg,1,1,ntmax), repmat(reshape(1:ntmax,1,1,[]),size(Ddist_gg,1),size(Ddist_gg,1),1), ampli.cov.parm(3), ampli.cov.parm(4), ampli.cov.parm(5), ampli.cov.parm(6), ampli.cov.parm(7));


% 4. Initizialization of the kriging weights and variance error
NEIGHG = nan(numel(path),neighg_nb);
LAMBDAG = NEIGHG;
NEIGHR = nan(numel(path),neighr_nb);
LAMBDAR = NEIGHR;
S = nan(numel(path),1);
% LATLONAT_i=[g.lat3D(path) g.lon3D(path) datenum(g.time3D(path))];
LL_i=sub2ind([g.nlat g.nlon],LAT(path),LON(path));
TT_i=AT(path);

% 5 Loop of scale for multi-grid path
for i_pt = 1:numel(path)
    
    % 1. Neighbor from grid
    % maybe faster?
    % neigh=find(Dtime_gg(TT_i(1:i_pt-1),TT_i(i_pt))<ampli.cov.parm(4)*wradius & Ddist_gg(LL_i(1:i_pt-1),LL_i(i_pt))<ampli.cov.parm(3)*wradius);

    neighg=find(Dtime_gg(TT_i(1:i_pt-1),TT_i(i_pt))<ampli.cov.parm(4)*wradius);
    neighg=neighg(Ddist_gg(LL_i(neighg),LL_i(i_pt))<ampli.cov.parm(3)*wradius);

    Cgp = ampli.cov.parm(2).*Gneiting(Ddist_gg(LL_i(neighg),LL_i(i_pt)), Dtime_gg(TT_i(neighg),TT_i(i_pt)), ampli.cov.parm(3), ampli.cov.parm(4), ampli.cov.parm(5), ampli.cov.parm(6), ampli.cov.parm(7));
    [Cgp,tmp]=maxk(Cgp,neighg_nb);
    neighg=neighg(tmp);
    
    Cgg = ampli.cov.parm(2).*Gneiting(Ddist_gg(LL_i(neighg),LL_i(neighg)), Dtime_gg(TT_i(neighg),TT_i(i_pt)), ampli.cov.parm(3), ampli.cov.parm(4), ampli.cov.parm(5), ampli.cov.parm(6), ampli.cov.parm(7));
    
    neighr=find(Ddist_gr(:,LL_i(i_pt))<ampli.cov.parm(3)*wradius & Dtime_gr(:,TT_i(i_pt))<ampli.cov.parm(4)*wradius);
    Crp = ampli.cov.parm(2).*Gneiting(Ddist_gr(neighr,LL_i(i_pt)), Dtime_gr(neighr,TT_i(i_pt)), ampli.cov.parm(3), ampli.cov.parm(4), ampli.cov.parm(5), ampli.cov.parm(6), ampli.cov.parm(7));
    [Crp,tmp]=maxk(Crp,neighr_nb);
    neighr=neighr(tmp);
    
    Crg = ampli.cov.parm(2).*Gneiting(Ddist_gr(neighr,LL_i(neighg)), Dtime_gr(neighr,TT_i(neighg)), ampli.cov.parm(3), ampli.cov.parm(4), ampli.cov.parm(5), ampli.cov.parm(6), ampli.cov.parm(7));
    Crr = ampli.cov.C(neighr,neighr);
    
    l =  [Cgg Crg'; Crg Crr] \  [Cgp; Crp];

    NEIGHG(i_pt,1:numel(neighg))=neighg;
    NEIGHR(i_pt,1:numel(neighr))=neighr;
    LAMBDAG(i_pt,1:numel(neighg)) = l(1:numel(neighg));
    LAMBDAR(i_pt,1:numel(neighr)) = l(numel(neighg)+1:end);
    
    S(i_pt) = ampli.cov.parm(1) + ampli.cov.parm(2) - l'*[Cgp; Crp];

end










Gneiting = @(dist,time,range_dist,range_time,delta,gamma,beta) 1./( (time./range_time).^(2.*delta) +1 ) .* exp(-(dist./range_dist).^(2.*gamma)./((time./range_time).^(2.*delta) +1).^(beta.*gamma) );


for i_l=1:numel(latlonmask)
    %gampli.An_est_tmp=nan(1,g.nat); gampli.An_sig_tmp=gampli.An_est_tmp;
    for i_t=1:g.nat
        neighr=gampli.Ddist_sf(:,latlonmask(i_l))<ampli.cov.parm(3)*2 & gampli.Dtime_sf(:,i_t)<ampli.cov.parm(4)*2;
        Cab = ampli.cov.parm(2).*Gneiting(gampli.Ddist_sf(neighr,latlonmask(i_l)), gampli.Dtime_sf(neighr,i_t), ampli.cov.parm(3), ampli.cov.parm(4), ampli.cov.parm(5), ampli.cov.parm(6), ampli.cov.parm(7));
        lambda = ampli.cov.C(neighr,neighr)  \  Cab;
        gampli.An_est(i_t,latlonmask(i_l)) = lambda' * ampli.An(neighr);
        gampli.An_sig(i_t,latlonmask(i_l)) = ampli.cov.parm(1) + ampli.cov.parm(2) - lambda' * Cab;
        %gampli.An_est_tmp(i_t) =  lambda' * gampli.An(id2);
        %gampli.An_sig_tmp(i_t) = ampli.cov.parm(1) + ampli.cov.parm(2) - lambda' * Cab;
    end
    %gampli.An_est(:,ll(i_l)) = gampli.An_est_tmp;
    %gampli.An_sig(:,ll(i_l)) = gampli.An_sig_tmp;
end

gampli.An_est=reshape(gampli.An_est',g.nlat,g.nlon,g.nat);
gampli.An_sig=reshape(gampli.An_sig',g.nlat,g.nlon,g.nat);