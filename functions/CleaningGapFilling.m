function d = CleaningGapFilling(d)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOT MAINTAINED !!! 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i_d=1:numel(d)
    
    % Replace to the lowest altitudinal layer of 0 with the first none-zero layer
    lev = find(all(isnan(d(i_d).dens))==0);
    for l=1:lev(1)-1
        d(i_d).dens(:,l)=d(i_d).dens(:,lev(1));
    end
    
    % replace lower layer if upper layer as value
    
    
    % replace close inside values
    
    % replace the NaN by 0.
    d(i_d).dens(isnan(d(i_d).dens))=0;
    
    % remove outliar
    % if any(strcmp(d(i_d).name,outliar))
    %     d2log=d(i_d).dens(:);
    %     d2log(d2log==0 | isnan(d2log))=[]; % Removing two 0 values...
    %     logpd = fitdist((d2log),'lognormal');
    %     prob = pdf(logpd,d(i_d).dens);
    %     id=prob>0&prob<.0001;
    %     d(i_d).dens(id)=0;
    % %   imagesc(datenum(d(i_d).time), d(i_d).interval*(1/2:d(i_d).maxrange),log(d(i_d).dens)','AlphaData',~id')
    % %   set(gca, 'YDir', 'normal'); xlabel('date'); ylabel('[m]');
    % end
    %
    
end





% for i_d=1:numel(d)
% t= flipud(log(d(i_d).dens'+0.01));
% imagesc(t)
%     t = (t-min(t(:)))/range(t(:))
%     rgbImage = ind2rgb(+1,colormap );
%     imshow(rgbImage);
%     imwrite(flipud(d(i_d).dens'),['figure/ToCorrect/' num2str(i_d) '_' d(i_d).name '.png'])
% end
%
% keyboard;



% Photoshop
% psnewdocmatrix(flipud(d(i_d).dens'));
% d(i_d).dens = flipud(psgetpixels()');




% for i_d=1:numel(d)
%     xg1=( 1:size(d(i_d).dens,1) )*10000;
%     xg2=( 1:size(d(i_d).dens,2) )*1;
%     [X1,X2]=ndgrid(xg1,xg2);
%     id=isnan(d(i_d).dens);
%     F = scatteredInterpolant(X1(~id),X2(~id),d(i_d).dens(~id));
%     
%     d(i_d).dens(id) = max(F(X1(id),X2((id))), 0);
%     
% end

% for i_d=1:nueml(d)
%     figure(1); clf; set(gcf, 'Color', 'w');
%     subplot(2,1,1); hold on;
%     imagesc(datenum(dini(i_d).time), dini(i_d).interval*(1/2:dini(i_d).levels), log(dini(i_d).dens)','AlphaData',~isnan(dini(i_d).dens'))
%     plot([datenum(start_date) datenum(end_date)],[dini(i_d).height dini(i_d).height],'r','linewidth',2); caxis([-5 5])
%     xlabel('date'); ylabel('[m]'); title('Density [bird/m^3]'); colorbar;
%     datetick('x');  set(gca, 'YDir', 'normal'); datetick('x');
%     set(gca,'XTick',datenum(dini(i_d).time(1):4/24:dini(i_d).time(end)))
%     set(gca,'XTickLabel',datestr(dini(i_d).time(1):4/24:dini(i_d).time(end),'dd-mm HH'))
%     axis([datenum(start_date) datenum(end_date) 0 5000]);
    
%     subplot(2,1,2); hold on;
%     i_dd=strcmp(dini(i_d).name,{d.name});
%     if any(i_dd)
%         imagesc(datenum(d(i_dd).time), d(i_dd).interval*(1/2:d(i_dd).levels), log(d(i_dd).dens)')
%         rectangle('Position',[datenum(start_date) 0 datenum(end_date)-datenum(start_date) d(i_dd).height],'FaceColor',[1 1 1])
%         plot([datenum(start_date) datenum(end_date)],[d(i_dd).height d(i_dd).height],'r','linewidth',2); caxis([-5 5])
%         xlabel('date'); ylabel('[m]'); title('Density [bird/m^3]'); colorbar;
%         datetick('x'); axis([datenum(start_date) datenum(end_date) 0 5000]); set(gca, 'YDir', 'normal')
%         
%     end
    
%     subplot(2,1,2); hold on;
%     imagesc(datenum(dini(i_d).time), dini(i_d).interval*(1/2:dini(i_d).levels), dini(i_d).DBZH','AlphaData',~isnan(dini(i_d).DBZH'))
%     plot(datenum(dini(i_d).time), dini(i_d).tcrw*5000/max( dini(i_d).tcrw),'r-')
%     plot(datenum(dini(i_d).time([1 end])),[1 1].*0.005*5000/max( dini(i_d).tcrw),'-k')
%     imagesc(datenum(dini(i_d).time), dini(i_d).interval*(1/2:dini(i_d).maxrange), -50*ones(size(dini(i_d).DBZH')),'AlphaData', 0.6*repmat(dini(i_d).tcrw'>0.005,size(dini(i_d).DBZH,2),1))
%     xlabel('date'); ylabel('[m]'); title('Reflectivity []'); colorbar;
%     datetick('x');  set(gca, 'YDir', 'normal')
%     set(gca,'XTick',datenum(dini(i_d).time(1):4/24:dini(i_d).time(end)))
%     set(gca,'XTickLabel',datestr(dini(i_d).time(1):4/24:dini(i_d).time(end),'dd-mm HH'))
%     axis([datenum(start_date) datenum(end_date) 0 5000]);
%     
%     % export_fig(['figure/corrected/' num2str(i_d) '_' d(i_dd).name '.png'],'-png')
% end

% MPS-NOT USED
% s = system('bash -c
% "/mnt/c/Users/rafnu/Documents/GitHub/G2S/build/c++-build/server"');
% serverAddress='localhost';
% path=getenv('PATH');
% newpath=strcat(path,'C:\Program Files\ZeroMQ 4.0.4\bin;');
% setenv('PATH',newpath);
% serverAddress='tesla-k20c.gaia.unil.ch'; % 'tesla-k20c.gaia.unil.ch'
% for i_d=1:numel(d)
%     a{i_d} = d(i_d).dens(:,1:20)'; 
% end
% for i_d=1:numel(d)
%     parm={'-sa',serverAddress,'-a','qs','-ti',a{i_d},'-di',a{i_d},'-dt',zeros(1,1),'-k',1,'-n',9,'-s',100,'-j',1};
%     data=g2s(parm{:});
% end
% subplot(2,1,1)
% imagesc(datenum(d(i_d).time), (1:20)*200, log(a),'AlphaData',~isnan(a))
% subplot(2,1,2)
% imagesc(datenum(d(i_d).time), (1:20)*200, log(data))



end