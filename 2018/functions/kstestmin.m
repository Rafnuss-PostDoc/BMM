function p = kstestmin(data,f,parm)
tmp = f(data,parm);
[~,p]=kstest((tmp-nanmean(tmp))./nanstd(tmp));
p=-p;
end