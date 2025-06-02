function s = computeShares(p,params)
    % Compute plant shares given vector or matrix of prices
    if params.formulation=="D"
        share_exp = exp(params.alpha*p.*(1-params.sigma));
    elseif params.formulation=="M"
        share_exp = exp(params.alpha*(p.*(1-params.sigma)-params.T));
    end
    s = share_exp./sum(share_exp);
end