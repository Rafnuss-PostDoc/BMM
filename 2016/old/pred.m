function out = pred(x,s,nb_poly,scoret,dateradar)


% p=x(1:nb_poly+1);
% ms=x(nb_poly+1+(1:numel(dc)));
% mtslope = x(nb_poly+1+(1:numel(dc)));
% 
% pred = @(x) polyval(p,scoret) .* radar*ms; 
% pred = @(x) polyval(p,scoret) .* radar*ms .* (radar*mtslope.*scoret+1);


p=x(1:nb_poly+1);
A = reshape(x(nb_poly+2:end),s(1),s(2));
out = polyval(p,scoret) + A(dateradar);

end