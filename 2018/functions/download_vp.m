function d = download_vp(start_date,end_date, quantity, folder_root)

% Find available data
websave('data/coverage.csv','https://lw-enram.s3-eu-west-1.amazonaws.com/coverage.csv');
T = readtable('data/coverage.csv');

% Create a cell for all data. 
cr = unique(T.countryradar);
clear d
for i_d = 1:numel(cr)
    d(i_d).name = cr{i_d};
    d(i_d).date = unique(char(T(strcmp(T.countryradar,d(i_d).name),:).date,'uuuu-MM'),'rows');
    
    tmp = datetime(d(i_d).date,'InputFormat','yyyy-MM');
    d(i_d).date(tmp<start_date | tmp>end_date,:)=[];
end

for i_d=1:numel(d)
    
    d(i_d).time = [];
%     d(i_d).stime = [];
%     d(i_d).etime = [];
    for i_quantity = 1:numel(quantity)
        d(i_d).(quantity{i_quantity}) = [];
    end
    
    d(i_d).levels=nan(size(d(i_d).date,1),1);

    for i_dd = 1:size(d(i_d).date,1)
        
        filename = [d(i_d).name, strrep(d(i_d).date(i_dd,:),'-','')];
        
        % Download all the available file and unzip them
        if exist([folder_root filename '.zip'], 'file') ~= 2
            path = ['https://lw-enram.s3-eu-west-1.amazonaws.com/', d(i_d).name(1:2), '/', d(i_d).name(3:5), '/', extractBefore(d(i_d).date(i_dd,:),'-'), '/' ];
            try
                websave([folder_root filename  '.zip'],[path filename '.zip']);
                unzip([folder_root filename  '.zip'],[folder_root filename '/']);
            catch
                delete([folder_root filename  '.zip'])
                disp(['Error in downloading or unziping : ' filename ' (' i_dd ')'])
            end
        else
            unzip([folder_root filename  '.zip'],[folder_root filename '/']);
        end

        % Get all the individual h5 files 
        files = dir([folder_root, filename, '\**\*.h5']);
        
        % Filter date out of range
        if numel(files)>0
            tmp = cellfun(@(x) datetime( x(10:17) ,'InputFormat','yyyyMMdd'),{files.name},'UniformOutput',0);
            files = files([tmp{:}]>start_date & [tmp{:}]<end_date );
        end
        
        % Get the information about the radar. Should be the same for all
        % files... to be checked...
        if i_dd==1 && numel(files)>0
            d(i_d).lat = h5readatt([files(1).folder '/' files(1).name ],'/where','lat');
            d(i_d).lon = h5readatt([files(1).folder '/' files(1).name ],'/where','lon');
            d(i_d).height = h5readatt([files(1).folder '/' files(1).name ],'/where','height');
            % CHANGING ! d(i_d).levels = h5readatt([files(1).folder '/' files(1).name ],'/where','levels');
            d(i_d).interval = h5readatt([files(1).folder '/' files(1).name ],'/where','interval');
            d(i_d).maxrange = h5readatt([files(1).folder '/' files(1).name ],'/how','maxrange');
            
            quantityMatched = mapQuantity(files(1),quantity);
        end
        
        if numel(files)>0
            d(i_d).levels(i_dd) = h5readatt([files(1).folder '/' files(1).name ],'/where','levels');
        else
            d(i_d).levels(i_dd) = 0;
        end
        
        q = repmat({nan(numel(files), d(i_d).levels(i_dd))},numel(quantity),1);
        time = NaT(numel(files), 1);
        stime = NaT(numel(files), 1);
        etime = NaT(numel(files), 1);
        
        for i_file = 1: numel(files)
            
            % Check if file empty
            try
                h5info([files(i_file).folder '/' files(i_file).name],'/dataset1/');
            catch
                continue
            end
            
            % clc;h5disp([files(i_file).folder '/' files(i_file).name])
            
%             % Time
%             date =  h5readatt([files(i_file).folder '/' files(i_file).name ],'/what','date');
%             starttime = h5readatt([files(i_file).folder '/' files(i_file).name ],'/how','starttime');
%             endtime = h5readatt([files(i_file).folder '/' files(i_file).name ],'/how','endtime');
%             
%             if str2double(starttime(end-2:end-1))>=60
%                starttime(end-4:end-3) = sprintf('%02i',str2double(starttime(end-4:end-3))+1);
%                starttime(end-2:end-1) = sprintf('%02i',mod(str2double(starttime(end-2:end-1)),60));
%             end
%             if str2double(endtime(end-2:end-1))>=60
%                endtime(end-4:end-3) = sprintf('%02i',str2double(endtime(end-4:end-3))+1);
%                endtime(end-2:end-1) = sprintf('%02i',mod(str2double(endtime(end-2:end-1)),60));
%             end
%             stime(i_file) = datetime( [date(1:end-1) starttime(1:end-1)] ,'InputFormat','yyyyMMddHHmmss');
%             etime(i_file) = datetime( [date(1:end-1) endtime(1:end-1)] ,'InputFormat','yyyyMMddHHmmss');
            
            if strfind(files(i_file).name(10:22),'T')>0
                time(i_file) = datetime( files(i_file).name(10:22) ,'InputFormat','yyyyMMdd''T''HHmm');
            else
                time(i_file) = datetime( files(i_file).name(10:21) ,'InputFormat','yyyyMMddHHmm');
            end
            
            
            for i_quantity = 1:numel(quantity)
                % attr = h5readatt([files(i_file).folder '/' files(i_file).name ],[quantityMatched{i_quantity} '/what'],'quantity');
                % if strcmp( attr(1:end-1), quantity{i_quantity})
                try 
                   qq  = h5read([files(i_file).folder '/' files(i_file).name ], [quantityMatched{i_quantity} '/data']);
                catch
                    quantityMatched = mapQuantity(files(i_file),quantity);
                    qq = h5read([files(i_file).folder '/' files(i_file).name ], [quantityMatched{i_quantity} '/data']);
                end
                
                if numel(qq) == size(q{i_quantity},2)
                    q{i_quantity}(i_file,:)=qq;
                elseif numel(qq) > size(q{i_quantity},2)
                    dl = numel(qq)-size(q{i_quantity},2);
                    q{i_quantity} = [q{i_quantity} nan(size(q{i_quantity},1), dl)  ];
                    q{i_quantity}(i_file,:) = qq;
                else
                    q{i_quantity}(i_file,:)=[qq  nan(1,size(q{i_quantity},2)-numel(qq)) ];
                end
                
                
            end
            
        end
        
        d(i_d).time=[d(i_d).time; time];
%         d(i_d).stime=[d(i_d).stime; stime];
%         d(i_d).etime=[d(i_d).etime; etime];
        
        for i_quantity = 1:numel(quantity)
            if size(q{i_quantity},2) == size(d(i_d).(quantity{i_quantity}),2)
                d(i_d).(quantity{i_quantity}) = [d(i_d).(quantity{i_quantity}) ; q{i_quantity} ];
            elseif size(q{i_quantity},2) > size(d(i_d).(quantity{i_quantity}),2)
                dl = size(q{i_quantity},2) - size(d(i_d).(quantity{i_quantity}),2);
                d(i_d).(quantity{i_quantity}) = [d(i_d).(quantity{i_quantity}) nan(size(d(i_d).(quantity{i_quantity}) ,1), dl)  ];
                d(i_d).(quantity{i_quantity}) = [d(i_d).(quantity{i_quantity}) ; q{i_quantity} ];
            else
                dl = size(d(i_d).(quantity{i_quantity}),2) - size(q{i_quantity},2);
                q{i_quantity} = [q{i_quantity} nan(size(q{i_quantity},1), dl)  ];
                d(i_d).(quantity{i_quantity}) = [d(i_d).(quantity{i_quantity}) ; q{i_quantity} ];
            end
        end
        
        % d(i_d).dens=append(d(i_d).dens, timeseries(dens,datestr(time)) );
        if exist([folder_root filename], 'dir') == 7
            rmdir([folder_root filename], 's')
        end
    end
end


% remove date
d = rmfield(d,'date');


% Remove nan values and -1000 values
idd=[];
for i_d=1:numel(d)
    
    for i_quantity=1:numel(quantity)
        d(i_d).(quantity{i_quantity})(d(i_d).(quantity{i_quantity})<-998) =NaN;
    end
    
    if numel(d(i_d).dens)==0
        idd=[idd;i_d];
    end
end
d(idd)=[];

end


function quantityMatched = mapQuantity(file,quantity)

% dataset string corresponding to the quantity desired.
quantityMatched = cell(size(quantity));

info = h5info([file.folder '/' file.name],'/dataset1/');
dataset1={info.Groups.Name};
for i_dataset = 1:numel(dataset1)
    attr = h5readatt([file.folder '/' file.name ],[dataset1{i_dataset} '/what'],'quantity');
    id =strcmp(quantity,attr(1:end-1));
    if any(id==1)
        quantityMatched{id} = dataset1{i_dataset};
    end
end
end