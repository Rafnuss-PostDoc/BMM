function quantityMatched = mapQuantity(file,quantity)

% dataset string corresponding to the quantity desired.
quantityMatched = cell(size(quantity));

info = h5info([file.folder '/' file.name],'/dataset1/');
dataset1={info.Groups.Name};
for i_dataset = 1:numel(dataset1)
    if contains(dataset1{i_dataset},'what')
        continue
    end
    attr = h5readatt([file.folder '/' file.name ],[dataset1{i_dataset} '/what'],'quantity');
    id = strcmp(quantity,attr(1:end-1));
    if any(id==1)
        quantityMatched{id} = dataset1{i_dataset};
    end
end
end