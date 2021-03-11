
% Located in 
function [e,h,M,S,A] = histvdens(X,W,binedge,norm)
    if ~exist('binedge','var'), binedge = 500; end
    if ~exist('norm','var'), norm = false; end
    id = ~isnan(X(:))&~isnan(W(:));
    X=X(id); W = W(id);
    if sum(id(:))==0
        e=nan;
        h=nan;
        M=nan;
        S=nan;   
        A=nan;
    else
        % [~,e,Y] = histcounts(X,edges);
        [Y,edge] = discretize(X,binedge,'IncludedEdge','right'); 
        [G,ID] = findgroups(Y);
        h=zeros(numel(edge),1);
        h(ID+1) = splitapply(@sum, W, G);  
        h = h(2:end); % h(1) is for all X outside the bins
        e = 0.5 * (edge(1:end-1) + edge(2:end));
        if norm, h=h./nansum(h); end
        M = X'*W/sum(W);
        S = sqrt(var(X,W));
        A = mean(W.*exp(1i*deg2rad(X)));
    end
end