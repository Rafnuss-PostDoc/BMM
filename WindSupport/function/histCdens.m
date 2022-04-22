function [b,h,M,S] = histCdens(X,W,binedge,norm)
    if ~exist('binedge','var'), binedge = -pi:0.1:pi; end
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
        [Y,e] = discretize(X,binedge,'IncludedEdge','right'); 
        [G,ID] = findgroups(Y);
        h=zeros(size(e));
        h(ID+1) = splitapply(@sum, W, G);  %W.*abs(X)
        h = h(2:end); % h(1) is for all X outside the bins
        b = 0.5 * (e(1:end-1) + e(2:end));
        if norm, h=h./nansum(h); end
        M = X'*W/sum(W);
        S = sqrt(var(X,W));
        %A = mean(W.*exp(1i*deg2rad(X)));
    end
end