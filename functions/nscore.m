function Nscore = nscore(z) % X,support_dist,method,extrapolationMethod,plotit

[f,xi]=ksdensity(z,'NumPoints',100);

prior = f./sum(f);
prior(prior<eps)=2*eps;
cdf = cumsum(prior) ./sum(prior);

Nscore.T_F = griddedInterpolant(xi,cdf,'pchip','pchip');
Nscore.Tinv_F = griddedInterpolant(cdf,xi,'pchip','pchip');

z_T_mean=nanmean(norminv(Nscore.T_F(z(:))));
z_T_std=nanstd(norminv(Nscore.T_F(z(:))));

%Nscore.inverse = @(y) Nscore.Tinv_F(normcdf(y)); % back-transform a value in normal space by taking the normcdf.
%Nscore.forward = @(y) norminv( Nscore.T_F(y) );
Nscore.inverse = @(y) reshape(Nscore.Tinv_F(normcdf(y(:)*z_T_std+z_T_mean)),size(y)); % back-transform a value in normal space by taking the normcdf.
Nscore.forward = @(y) reshape((norminv(Nscore.T_F(y(:)))-z_T_mean)/z_T_std,size(y));

% kernel_y_ns = norminv(cdf);
% kernel_y_ns_mid = ( kernel_y_ns(1:end-1)+kernel_y_ns(2:end) ) /2;
% Nscore.dist = @(mu,sigma) [ normcdf(kernel_y_ns_mid,mu,sigma) ; 1] - [0 ; normcdf(kernel_y_ns_mid(),mu,sigma)]; 


end