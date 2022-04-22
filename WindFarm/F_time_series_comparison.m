ts=readtable('WT_data\time_series_60min_singleindex.csv');
ts.datetime = datetime(ts.cet_cest_timestamp,'InputFormat','yyyy-MM-dd''T''HH:mm:SSZ','TimeZone', 'local');
ts(year(ts.datetime)~=2018,:)=[];

tmp = energy;
tmp(repmat(~inCountry(:,:,2),1,1,numel(time)))=0;
time.TimeZone='local';

figure; hold on
plot(ts.datetime,ts.DE_wind_onshore_generation_actual*1e6*60*60) % MW ->J
plot(time,reshape(nansum(nansum(tmp),2),1,[]))
legend('DE_wind_onshore_generation_actual','estimated from our model')

figure;
plot(ts.datetime, ts.DE_wind_onshore_generation_actual*1e6*60*60./reshape(nansum(nansum(tmp),2),1,[])')
ylabel('Efficiency=ratio of actual power / power potential')

figure; hold on
plot(hour(ts.datetime ),ts.DE_LU_price_day_ahead,'.k')
plot(0:23,splitapply(@nanmean,ts.DE_LU_price_day_ahead,hour(ts.datetime )+1),'or')
ylabel('Price in euro'); xlabel('hour of day')

figure; hold on;
plot(time,reshape(nansum(nansum(birdAtRisk,1),2),1,[]))
yyaxis 'right'
plot(ts.datetime,ts.DE_LU_price_day_ahead)

figure;
plot(day(ts.datetime,'dayofyear')+hour(ts.datetime)/24,ts.DE_LU_price_day_ahead)