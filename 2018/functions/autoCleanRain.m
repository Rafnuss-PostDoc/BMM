function rain = autoCleanRain(di,rank1,thr1,rank2,thr2,win2)

% Remove ground scattering
di.eta(:,1:di.scatter_lim) = nan; 
di.DBZH(:,1:di.scatter_lim_DBZH) = nan; 

% deoft has a lot of noise in eta. Take dens and clean it. 
% i_d=strcmp({dc.name},'depro'); %deoft, dehnr, depro
% di.dens2=di.dens;
% di.dens2(all(isnan(di.dens3),2),:) = NaN;
% for i_t=1:numel(di.time)
%     firt=find(di.dens(i_t,:)>0,1,'last');
%     id=di.dens2(i_t,1:firt)==0 | isnan(di.dens2(i_t,1:firt));
%     di.dens2(i_t,id)= di.dens3(i_t,id);
% end

rain = false(size(di.DBZH,1),1);

% DBZH>-1 db for the 8th highest value 
a=sort(di.DBZH,2,'descend','MissingPlacement','last');
rain(a(:,rank1)>thr1)=true;

% add padding of 30min
% rain = movmax(di.rain,13);

% Hytherisis: if DBZH>-6 db for the 5th highest value and within 13*5min to
% previous rain
id = 1;
while sum(id)>0
   id = a(:,rank2)>thr2 & movmean(rain, win2)>0 & ~rain;
   rain(id)=true;
end

% Remove for data which last 30minutes or less
n=30/5;
B1 = ones(2*n+1,1); B1(n+1)=0;
id_nan = di.day | rain | all(isnan(di.eta),2);
rain(~id_nan & conv2(id_nan,B1,'same') > n)=true;

end