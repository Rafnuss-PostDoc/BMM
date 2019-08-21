function out = pred2(x,n_dc,n_uniqueDate,nb_poly,data)


% p=x(1:nb_poly+1);
% ms=x(nb_poly+1+(1:numel(dc)));
% mtslope = x(nb_poly+1+(1:numel(dc)));
% 
% pred = @(x) polyval(p,scoret) .* radar*ms; 
% pred = @(x) polyval(p,scoret) .* radar*ms .* (radar*mtslope.*scoret+1);

S = reshape(x(1:n_dc),n_dc,1);
A = reshape(x(n_dc+(1:n_dc*n_uniqueDate)),n_dc,n_uniqueDate);
p = x(end-nb_poly:end);

out = S(data.i_d)+ A(data.dateradar) + polyval(p,data.scoret);

end