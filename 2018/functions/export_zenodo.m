load('data/dc_corr')
load('./data/Density_inference.mat'); 
load('./data/BelowRadarMPS.mat');

% Run the begining of `Density_inference.mlx` up to "% Remove time with all
% day or no data to save space" to generate `data` and `data_denss_sim` 
% before time without data is removed. 

%% Export Zenodo
dexp = dc;
kep=false(numel(dc),numel(dc(1).time));
for i_d=1:numel(dexp)
    dexp(i_d).dens=dc(i_d).dens4;
    dexp(i_d).dens(:,1:dc(i_d).scatter_lim-1)=nan;
    dexp(i_d).densSIMmean = mean(MPS{i_d},3);
    
    dexp(i_d).sd_vvp=dexp(i_d).sd_vvp2;
    dexp(i_d).volDir=dexp(i_d).VolBelow;
    
    dexp(i_d).uw = real(dexp(i_d).ws);
    dexp(i_d).vw = imag(dexp(i_d).ws);
    
    dexp(i_d).ubs = data.vs(:,i_d);
    dexp(i_d).vbs = data.us(:,i_d);
    dexp(i_d).denss = mean(data_denss_sim(:,i_d,:),3);
    
    kep(i_d,:) = any(~isnan(dexp(i_d).dens'));
end

dexp = rmfield(dexp,{'DBZH','ff','dd','interval','scatter_lim','levels','dusk','dawn','sunset','sunrise','day','sd_vvp2','eta','u','v','dens2','dens3','dens4','time','alt','maxrange','VolBelow','ws','wt'});

for i_d=1:numel(dexp)
    
    dexp(i_d).dens = round(dexp(i_d).dens(any(kep),:),2)';
    dexp(i_d).sd_vvp = round(dexp(i_d).sd_vvp(any(kep),:),2)';
    dexp(i_d).ub = round(dexp(i_d).ub(any(kep),:),2)';
    dexp(i_d).vb = round(dexp(i_d).vb(any(kep),:),2)';
    dexp(i_d).ui = round(dexp(i_d).ui(any(kep),:),2)';
    dexp(i_d).vi = round(dexp(i_d).vi(any(kep),:),2)';
    dexp(i_d).uw = round(dexp(i_d).uw(any(kep),:),2)';
    dexp(i_d).vw = round(dexp(i_d).vw(any(kep),:),2)';
    dexp(i_d).insect = round(dexp(i_d).insect(any(kep),:),2)';
    dexp(i_d).bird = round(dexp(i_d).bird(any(kep),:),2)';
    dexp(i_d).densSIMmean = round(dexp(i_d).densSIMmean(any(kep),:),2)';
    dexp(i_d).ubs = round(dexp(i_d).ubs(any(kep)),2)';
    dexp(i_d).vbs = round(dexp(i_d).vbs(any(kep)),2)';
    dexp(i_d).denss = round(dexp(i_d).denss(any(kep)),2)';
    
    fileID = fopen(['data/zenodo/dc_' dexp(i_d).name '.json'],'w');
    fprintf(fileID,jsonencode(dexp(i_d)));
    fclose(fileID);
end

fileID = fopen('data/zenodo/time.json','w');
fwrite(fileID,jsonencode(dc(1).time(any(kep))));
fclose(fileID)
