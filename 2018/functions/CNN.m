load('data/dc_corr.mat')

i_d=4;
dc2 = dc;


dc2(i_d).dens = log10(dc(i_d).eta/11);
dc2(i_d).dens(dc(i_d).day,:)=nan;
dc2(i_d).dens(:,1:dc(i_d).scatter_lim)=nan;

dc2(i_d).dens2 = log10(dc(i_d).dens2);
dc2(i_d).dens2(:,1:dc(i_d).scatter_lim)=nan;


dc2(i_d).cat = ones(size(dc2(i_d).sd_vvp));
% cat 3: value in eta different from dens2
dc2(i_d).cat((dc2(i_d).dens2-dc2(i_d).dens)>0.000001) = 3;
% cat 1: value in eta -> nan or 0 in dens2
dc2(i_d).cat(isnan(dc2(i_d).dens2) & ~isnan(dc2(i_d).dens)) = 2;




X = cat(3, dc2(i_d).dens, dc2(i_d).DBZH, dc2(i_d).sd_vvp);

y = dc2(i_d).cat;



dc2(i_d).dens3 = dwt2(dc2(i_d).dens,'db5');
















%%

i=1;

idx=([datasample(find(y==1),2000);datasample(find(y==2),2000);datasample(find(y==3),2000)]);
TrainX=nan(51,5,3,length(idx));
TrainY=nan(1,length(idx));
[p1Array,p2Array]=ind2sub(size(y),idx);

while i<=size(TrainX,4)
    
    p1=p1Array(i);
    p2=p2Array(i);
    try
        TrainX(:,:,:,i)=X(p1-25:p1+25,p2-2:p2+2,:);
        TrainY(i)=y(p1,p2);
        
    catch
       
    end
    i=i+1;
end
TrainX(:,:,:,isnan(TrainY))=[];
TrainY(isnan(TrainY))=[];

perm=randperm(length(TrainY));
TrainX=TrainX(:,:,:,perm);
TrainY=TrainY(perm);

%%
layers = [ ...
    imageInputLayer([51 5 3])
    convolution2dLayer(5,20)
    fullyConnectedLayer(50*5*3 * 3)
%     reluLayer
%     fullyConnectedLayer(50*5*3 * 10)
    reluLayer
    fullyConnectedLayer(50*5*3 * 1)
    reluLayer
    fullyConnectedLayer(200)
    fullyConnectedLayer(3)
    softmaxLayer
    classificationLayer];

options = trainingOptions('adam', ...
    'MaxEpochs',20,...
    'InitialLearnRate',1e-6, ...
    'Verbose',false, ...
    'Plots','training-progress');
%%
net = trainNetwork(TrainX,categorical(TrainY,[1 2 3]),layers,options);
