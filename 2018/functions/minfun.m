function r = minfun(x,d,tn,f_tid,bejab_rain)
th1=x(1);
buf1 = round(x(2));
q=x(3);
th2 =x(4);
buf2=round(x(5));
rain = false(size(d.time));

rain(nanmedian(d.DBZH,2)>th1)=true;

if buf1~=0
    rainCS = false(tn);
    id = f_tid(d.time);
    rainCS(id)=rain;
    rainCS = movmax(rainCS,buf1);
    rain = rainCS(id)';
end

a=quantile(d.sd_vvp,q,2)<th2;
rain(a)=true;

if buf2~=0
    rainCS = false(tn);
    id = f_tid(d.time);
    rainCS(id)=rain;
    rainCS = movmax(rainCS,buf2);
    rain = rainCS(id)';
end

r = [ sum(rain&bejab_rain) sum(~rain&bejab_rain); sum(rain&~bejab_rain) sum(~rain&~bejab_rain)];

rr=sum(r(2)*2+r(3));
