%% Export Zenodo
dexp = dc;
kep=false(numel(dc),numel(dc(1).time));
for i_d=1:numel(dexp)
    dexp(i_d).dens=dc(i_d).dens4;
    dexp(i_d).densSIMmean = mean(MPS{i_d},3);
    
    dexp(i_d).u=dexp(i_d).u2;
    dexp(i_d).v=dexp(i_d).v2;
    dexp(i_d).sd_vvp=dexp(i_d).sd_vvp2;
    dexp(i_d).volDir=dexp(i_d).VolBelow;
    
    kep(i_d,:) = any(~isnan(dexp(i_d).dens'));
end

dexp = rmfield(dexp,{'DBZH','ff','dd','interval','scatter_lim','levels','dusk','dawn','sunset','sunrise','day','sd_vvp2','eta','u2','v2','dens2','dens3','dens4','time','alt','maxrange','VolBelow'});

for i_d=1:numel(d_exp)
    
    dexp(i_d).dens = round(dexp(i_d).dens(any(kep),:),2)';
    dexp(i_d).sd_vvp = round(dexp(i_d).sd_vvp(any(kep),:),2)';
    dexp(i_d).u = round(dexp(i_d).u(any(kep),:),2)';
    dexp(i_d).v = round(dexp(i_d).v(any(kep),:),2)';
    dexp(i_d).windu = round(dexp(i_d).windu(any(kep),:),2)';
    dexp(i_d).windv = round(dexp(i_d).windv(any(kep),:),2)';
    dexp(i_d).insect = round(dexp(i_d).insect(any(kep),:),2)';
    dexp(i_d).densSIMmean = round(dexp(i_d).densSIMmean(any(kep),:),2)';
    
    fileID = fopen(['data/zenodo/dc_' dexp(i_d).name '.json'],'w');
    fprintf(fileID,jsonencode(dexp(i_d)));
    fclose(fileID);
end

fileID = fopen('data/zenodo/time.json','w');
fwrite(fileID,jsonencode(dc(1).time(any(kep))));
fclose(fileID)
