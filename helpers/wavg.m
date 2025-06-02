function p_avg = wavg(p,q)
% Function to compute weighted average of prices from county-level averages
    p_avg = sum(q.*p)/sum(q);
end