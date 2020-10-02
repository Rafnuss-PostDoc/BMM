function X=trend(t,lat,lon,poly_t_degree,poly_d_degree)
% assert that t, lat and lon have the same size

X=bsxfun(@power,t(:),poly_t_degree:-1:1);
for i1=poly_d_degree:-1:1
    X = [X lat(:).^i1];
end
for i2=poly_d_degree:-1:1
    X = [X lon(:).^i2];
end
X = [X ones(size(X,1),1)];
end