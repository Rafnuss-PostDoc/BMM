%% Export data for vpexplore https://bmm.raphaelnussbaumer.com/explorevp/


%% Export Clean 2018
load('data/dc_corr.mat'); load('./data/BelowRadarMPS.mat');

for i_d=1:numel(dc)
    dc(i_d).dens=dc(i_d).dens4;
    dc(i_d).dens(:,1:dc(i_d).scatter_lim-1) = mean(MPS{i_d}(:,1:dc(i_d).scatter_lim-1,:),3);
    
    dc(i_d).u=dc(i_d).u2;
    dc(i_d).v=dc(i_d).v2;
    dc(i_d).sd_vvp=dc(i_d).sd_vvp2;
end

d_zen = rmfield(dc,{'DBZH','ff','dd','interval','levels','dusk','dawn','sunset','sunrise','day','sd_vvp2','eta','u2','v2'});
d_exp = rmfield(d_zen,{'u','v','heightDEM','dens3','dens2','dens4','VolBelow','windu','windv','insect'});

for i_d=1:numel(d_exp)
    
    d_exp(i_d).time = datestr(d_exp(i_d).time,'yyyy-mm-dd HH:MM');
    d_exp(i_d).dens = round(d_exp(i_d).dens,2)';
    
    fileID = fopen(['data/vp/vp-clean/dc_' d_exp(i_d).name '.json'],'w');
    fprintf(fileID,jsonencode(d_exp(i_d)));
    fclose(fileID);
end



%% Export raw data 2018
load('data/d2018all.mat');
for i_d=1:numel(d)
    d(i_d).alt=d(i_d).interval*(1:d(i_d).levels)-d(i_d).interval/2;
end

d_raw = rmfield(d,{'interval','levels','eta','dd','ff','DBZH','sd_vvp'});

for i_d=1:numel(d_raw)
    d_raw(i_d).time = datestr(d_raw(i_d).time,'yyyy-mm-dd HH:MM');
    d_raw(i_d).dens = round(d_raw(i_d).dens,2)';
    %d_raw(i_d).ff = round(d_raw(i_d).ff,2);
    %d_raw(i_d).dd = round(d_raw(i_d).dd,2);
    %d_raw(i_d).DBZH = round(d_raw(i_d).DBZH,2);
    %d_raw(i_d).eta = round(d_raw(i_d).eta,2);
    %d_raw(i_d).sd_vvp = round(d_raw(i_d).sd_vvp,2);
    
    
    fileID = fopen(['data/vp/vp-raw/dc_' d_raw(i_d).name '.json'],'w');
    fprintf(fileID,jsonencode(d_raw(i_d)));
    fclose(fileID);
end

%% Export 2019
load('data/d2019all.mat');
for i_d=1:numel(d)
    d(i_d).alt=d(i_d).interval*(1:d(i_d).levels)-d(i_d).interval/2;
end

d_2019 = rmfield(d,{'interval','levels','eta','dd','ff','DBZH','sd_vvp'});

for i_d=1:numel(d_2019)
    d_2019(i_d).time = datestr(d_2019(i_d).time,'yyyy-mm-dd HH:MM');
    d_2019(i_d).dens = round(d_2019(i_d).dens,2)';
    %d_raw(i_d).ff = round(d_raw(i_d).ff,2);
    %d_raw(i_d).dd = round(d_raw(i_d).dd,2);
    %d_raw(i_d).DBZH = round(d_raw(i_d).DBZH,2);
    %d_raw(i_d).eta = round(d_raw(i_d).eta,2);
    %d_raw(i_d).sd_vvp = round(d_raw(i_d).sd_vvp,2);
    
    
    fileID = fopen(['data/vp/vp-2019/dc_' d_2019(i_d).name '.json'],'w');
    fprintf(fileID,jsonencode(d_2019(i_d)));
    fclose(fileID);
end

%% Radar name
radarname = unique({d_raw.name d_exp.name d_2019.name});
d_list=struct();
for i=1:numel(radarname)
    d_list(i).name = radarname{i};
        
    i_d=find(strcmp(d_list(i).name, {d_raw.name}));
    if numel(i_d)>0
        d_list(i).lat = d_raw(i_d).lat;
        d_list(i).lon = d_raw(i_d).lon;
        d_list(i).height = d_raw(i_d).height;
        d_list(i).maxrange = d_raw(i_d).maxrange;
        d_list(i).dt = minutes(median(diff(d_raw(i_d).time)));
        d_list(i).nanrawday2018 = sum(~all(isnan(d_raw(i_d).dens)|d_raw(i_d).dens==0,2))*5/60/24;
    end
    
    i_d=find(strcmp(d_list(i).name, {d_exp.name}));
    if numel(i_d)>0
        d_list(i).scatter_lim = d_exp(i_d).scatter_lim*200;
        d_list(i).nancleanday2018 = sum(~all(isnan(dc(i_d).dens3),2))*5/60/24;
    end
    
    i_d=find(strcmp(d_list(i).name, {d_2019.name}));
    if numel(i_d)>0
        d_list(i).lat = d(i_d).lat;
        d_list(i).lon = d(i_d).lon;
        d_list(i).height = d(i_d).height;
        d_list(i).maxrange = d(i_d).maxrange;
        d_list(i).dt = minutes(median(diff(d(i_d).time)));
        d_list(i).nanrawday2019 = sum(~all(isnan(d_2019(i_d).dens)|d_2019(i_d).dens==0,2))*5/60/24;
    end
end

fileID = fopen('data/vp/radar_list.json','w');
fprintf(fileID,jsonencode(d_list));
fclose(fileID);

%% Export Clean 2016
load('../2016/1-Cleaning/data/dc_corr.mat'); 


%% Figure 
figure('position',[0 0 800 600]);
h=worldmap([35 69], [-10 32]);  
setm(h,'frame','on','grid','off'); set(findall(h,'Tag','MLabel'),'visible','off'); set(findall(h,'Tag','PLabel'),'visible','off')
files={'gt30w020n90','gt30e020n90','gt30e020n40','gt30w020n40'};
for i=1:numel(files)
    [X,cmap,R] = geotiffread(['data/' files{i}]);
    geoshow(double(X),cmap,'DisplayType','texturemap')
end
demcmap(X)
scatterm([d_2019.lat], [d_2019.lon],100,'filled','MarkerFaceColor',[.898, .239, 0],'MarkerEdgeColor',[.82, .22, 0]);
scatterm([dc.lat],[dc.lon],100,'filled','MarkerFaceColor',[1,.914,0], 'MarkerEdgeColor',[.91,.831,0]);
scatterm([d_exp.lat], [d_exp.lon],100,'filled','MarkerFaceColor',[0.9255    0.4588    0.0196],'MarkerEdgeColor',[0.9255    0.4588    0.0196]);

print('CoverageMap3','-dpng','-r600')
