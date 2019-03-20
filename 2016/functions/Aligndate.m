function d_out = Aligndate(d,start_date,end_date,ddate,q)



d_out.time=start_date:ddate:end_date;
% Radar quantity
for i_q=1:numel(q.r)
    d_out.(q.r{i_q})=nan(numel(d),numel(d_out.time)-1);
    for i_date = 1:numel(d_out.time)-1
        for i_d=1:numel(d)
            id = d(i_d).time>=d_out.time(i_date) & d(i_d).time<d_out.time(i_date+1);
            weight=diff(double([0 max(0, (1:d(i_d).levels).*d(i_d).interval - d(i_d).height )]))/1000;
            d_out.(q.r{i_q})(i_d,i_date) =  nanmean(d(i_d).(q.r{i_q})(id,:),1)*weight';
        end
    end
end
% Altitude quantity
for i_q=1:numel(q.a)
    d_out.(q.a{i_q})=nan(numel(d),numel(d_out.time)-1,size(d(1).(q.a{i_q}),2));
    for i_date = 1:numel(d_out.time)-1
        for i_d=1:numel(d)
            id = d(i_d).time>=d_out.time(i_date) & d(i_d).time<d_out.time(i_date+1);
            d_out.(q.a{i_q})(i_d,i_date,:) =  nanmean(d(i_d).(q.a{i_q})(id,:));
        end
    end
end
% Surface quantity
for i_q=1:numel(q.s)
    d_out.(q.s{i_q})=nan(numel(d),numel(d_out.time)-1);
    for i_date = 1:numel(d_out.time)-1
        for i_d=1:numel(d)
            id = d(i_d).time>=d_out.time(i_date) & d(i_d).time<d_out.time(i_date+1);
            d_out.(q.s{i_q})(i_d,i_date) =  nanmean(d(i_d).(q.s{i_q})(id));
        end
    end
end


d_out.time=d_out.time(1:end-1);
d_out.lat =[d.lat]; 
d_out.lon=[d.lon]; 
d_out.maxrange = [d.maxrange];
d_out.n=numel(d_out.lat); 
d_out.height=[d.height];
d_out.name={d.name};
 

if 0
    figure; imagesc(1:numel(d),datenum(date15),d15');
    xlabel('Radars'); ylabel('Date')
    set(gca, 'YDir', 'normal');
    xticks(1:numel(d)); xticklabels({d.name}); xtickangle(90); datetick('y'); axis tight
    title('Density')
    
    figure;
    plot(datenum(date15),d15');
    legend({d.name})
    xlabel('Date'); ylabel('Density [bird/m^2]');
    set(gca, 'YDir', 'normal');
    datetick('x','dd-mm'); axis tight
    
    
    figure; load coastlines;
    subplot(1,2,1);worldmap('Europe'); plotm(coastlat, coastlon)
    scatterm(lat,lon,maxrange*4,nanmean(dhour'),'filled');
    title('Mean')
    subplot(1,2,2);worldmap('Europe'); plotm(coastlat, coastlon)
    scatterm(lat,lon,maxrange*4,nanstd(dhour'),'filled');
    title('std')
end

end