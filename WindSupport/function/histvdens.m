function [b,h,M,S] = histvdens(X,W,binedge,norm)
    if ~exist('binedge','var'), binedge = 0:0.1:30; end
    if ~exist('norm','var'), norm = false; end
    id = ~isnan(X(:))&~isnan(W(:));
    X=X(id); W = W(id);
    if sum(id(:))==0
        b=nan;
        h=nan;
        M=nan;
        S=nan; 
    else
        % [~,e,Y] = histcounts(X,edges);
        [Y,edge] = discretize(X,binedge,'IncludedEdge','right'); 
        [G,ID] = findgroups(Y);
        h=zeros(numel(edge),1);
        h(ID+1) = splitapply(@sum, W, G);  
        h = h(2:end); % h(1) is for all X outside the bins
        b = 0.5 * (edge(1:end-1) + edge(2:end));
        if norm, h=h./nansum(h); end
        M = X'*W/sum(W);
        S = sqrt(var(X,W));
    end
end