function out = pred3(x,n_dc,n_uniqueDate,nb_poly,data,dc)


% p=x(1:nb_poly+1);
% ms=x(nb_poly+1+(1:numel(dc)));
% mtslope = x(nb_poly+1+(1:numel(dc)));
% 
% pred = @(x) polyval(p,scoret) .* radar*ms; 
% pred = @(x) polyval(p,scoret) .* radar*ms .* (radar*mtslope.*scoret+1);

tfit = x(1:2);
t = [dc.lat]'*tfit(1)+tfit(2);
A = reshape(x(2+(1:n_dc*n_uniqueDate)),n_dc,n_uniqueDate);
p = x(end-nb_poly:end);

out = t(data.i_r)+ A(data.dateradar) + polyval(p,data.scoret);

end