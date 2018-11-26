function p = kstestmin(i,x)
tmp = x.^i;
[~,p]=kstest((tmp-nanmean(tmp))./nanstd(tmp));
p=-p;
end