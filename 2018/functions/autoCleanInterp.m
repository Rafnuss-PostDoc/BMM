function interp = autoCleanInterp(dens3,thr3,slo3)


% Fours elements are interpolated
% 1. Interpolate noise defined by mean of direct neigh < 1/2 elmt
B2 = 1/8*ones(3,3); B2(2,2)=0;
tmp = log10(dens3+1);
id_nan = isnan(tmp);
tmp2=conv2(~id_nan,B2,'same');
tmp(id_nan)=0;
id_noise = tmp > 2*conv2(tmp,B2,'same')./tmp2;

% 4. Rain on top of the sky which stop quickly
tmp = max(dens3(:,di.scatter_lim:di.scatter_lim+5),[],2)<dens3(:,(di.scatter_lim+5+1):end);
id_raintop = [false(size(dens3,1),di.scatter_lim+5) tmp];

% 2 and 3. Interpolate up to 5000m and the gaps
tmp = ~all(isnan(dens3),2);
dens3(tmp,25)=0;
id_top = isnan(dens3);
id_top(~tmp,:)=false;
id_top(:,1:25<=di.scatter_lim) = false;

% interp = id_noise|id_top|id_raintop;
interp = id_raintop;

end