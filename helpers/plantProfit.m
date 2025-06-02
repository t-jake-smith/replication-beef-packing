function Pi_plant = plantProfit(p,s,scen_results,params)
% COMPUTE PLANT LEVEL PROFITS

% Pull params and variables out of inputs
Q_n = params.Q_n;
X = scen_results.X;
X_n = sum(X);
M_n = Q_n - X_n;
b = params.b;
T = params.T;
pct_beef = params.pct_beef;

% County average price
p_n = sum(s.*p);

live_weight = 1250; % Assume avg cattle is 1250 lbs/head

% Spot market revenue (net of cattle and transportation costs)
Rev_spot = sum(M_n.*s.*(b-p-T),2);

% Contracted revenue (net of cattle and transportation costs)
% Compute p_bar: average price for each county
p_bar = zeros(size(Q_n));
for n=1:length(p_bar)
    p_bar_counties = params.reporting_regions(n,:);
    M_r = sum(M_n(p_bar_counties));
    p_bar(n) = sum((M_n(p_bar_counties)./M_r).*p_n(p_bar_counties),2);
end
Rev_ama = sum(X.*(b-p_bar-T),2);

% Total processing cost - intergrate over marginal cost
Q_plant = sum(s.*M_n + X,2);
K_plant = params.K;
C = zeros(size(Q_plant));
% Loop over plants to do integration
for i=1:length(C)
    q = linspace(0,Q_plant(i),1000);
    C(i) = trapz(q,params.marginal_cost(q,K_plant(i)));
end

% Profit
Pi_plant = (Rev_spot + Rev_ama - C)*pct_beef*live_weight;
end
