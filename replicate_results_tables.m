%% Parameters and Setup

% Add subfolders to the search path
addpath('helpers\');
addpath('data\');

% Set up struct for model parameters passed to functions
params = struct();

% Model parameters
alpha = 56.28;          % Feedlot price sensitivity parameter
c0 = 0.3608;            % Base marginal cost
pct_beef = 0.417;       % 1 lb live cattle = 0.417 lb beef
d_max = 650;            % Maximum distance plants bid on cattle
operation_days = 240;   % Full operating days per year for plants
supply_delta = 0.15;    % Change in supply for supply increase and decrease scenarios
e_d = 0.5;              % Elasticity of demand for beef
b = 3.977;              % Beef price ($/lb beef)
t = 1.6923e-4/pct_beef; % Transportation cost ($/(lb beef x mile))

% Economic parameters
params.alpha = alpha*pct_beef;  % Convert alpha to $/lb live cattle because that's what feedlot operators will respond to
params.c0 = c0;                 % Base marginal cost
params.lambda = 1.5;            % Marginal cost curvature above K
params.gamma = (0.5*c0)./0.2^params.lambda; % Marginal cost coefficient above K - calibrated so that 20% above capacity results in 50% higher MC
params.pct_beef = pct_beef;     % Need to save pct_beef to compute profits

%% Load data tables

% County-level fed cattle production data and county centroid locations
tbl_county = readtable("fed_cattle_data.csv");

% Plant locations and capacities
tbl_plant = readtable("plant_data.csv");
% Multiply plant capacity (head/day) by opearating days per year to get
% annual capacity
tbl_plant.K = tbl_plant.capacity*operation_days;

%% Define scenarios and import raw algorithm output for each scenario
% Scenario 1: Base model calibrated to 2022 data
% Scenario 2: Plant entry: New 2000 head/day plant in Glenwood, IA
% Scenario 3: Plant exit: Shutdown of 6000 head/day plant in Holcomb, KS
% Scenario 4: Supply increase: Exogenous cattle supply increase of 10%,
%             beef price decrease by 10% (due to move along demand curve)
% Scenario 5: Supply decrease: Exogenous cattle supply decrease of 10%,
%             beef price increase by 10% (due to move along demand curve)
% Scenario 6: Independent plants
% Scenario 7: Limit contracting: Max 50% contracting share in each region
% Scenario 8: Mill price
% Scenario 9: Collusion among top 4 firms

% Define scenario names for reporting results
scen_name = {"base","entry","exit","supply_inc","supply_dec","ind_plants","limit_amas","mill","collusion"};

% Initialize struct for saving results
results = struct();
results.T_scenarios = table(); % Initalize table for comparison of baseline and counterfactual policy scenarios

for scenario = 1:9
    %% Define scenario parameters
    % Start with baseline values
    scen_tbl_county = tbl_county;
    scen_tbl_plant = tbl_plant;
    scen_b = b;
    scen_formulation = "D";
    
    % Adjust parameters based on the scenario
    name = scen_name{scenario};

    % Plant entry
    if scenario==2
        % Add table row for entrant plant
        new_plant = {50,26,2000,41.153716384192684,-95.83194475181442,2000*operation_days};
        scen_tbl_plant = [tbl_plant;new_plant];
    
    % Plant exit
    elseif scenario==3
        % Remove table row for exiting plant
        scen_tbl_plant = tbl_plant([1:12,14:end],:);
    
    % Supply increase
    elseif scenario==4
        scen_tbl_county.fed_cattle = tbl_county.fed_cattle*(1+supply_delta);
        scen_b = b*(1-supply_delta/e_d);
    
    % Supply decrease
    elseif scenario==5
        scen_tbl_county.fed_cattle = tbl_county.fed_cattle*(1-supply_delta);
        scen_b = b*(1+supply_delta/e_d);
    
    % Independent plants
    elseif scenario==6
       scen_tbl_plant.firm_id = (1:height(scen_tbl_plant))';

    % Limit contracting
    elseif scenario==7
        scen_tbl_county.contract_share = min(tbl_county.contract_share,0.5);

    % Mill price formulation
    elseif scenario==8
        scen_formulation = "M";
    
    % Collusion among top 4 firms: change firm IDs so they act as one firm
    elseif scenario==9
        scen_tbl_plant.firm_id_save = tbl_plant.firm_id;        % Save original firm IDs
        scen_tbl_plant.firm_id(scen_tbl_plant.firm_id<=4) = 1;  % Change firm IDs for top 4 firms to 1
    end

    % Add beef price and formulation for scenario to params
    params.b = scen_b;                      % Beef price
    params.formulation = scen_formulation;  % Formulation (price discrimination or mill price)
    
    % Add data (or data augmented for counterfactual scenario) to params
    params.Q_n = (scen_tbl_county.fed_cattle)';     % County-level data
    params.K = scen_tbl_plant.K;                    % Plant capacities
    
    % Variables computed from data
    params.N = height(scen_tbl_county);             % Number of counties
    params.J = height(scen_tbl_plant);              % Number of plants
    params.D = distances((scen_tbl_county.lat)', ...% JxN distance matrix D, where D(j,n)=distance from j to n in miles
        (scen_tbl_county.lon)', ...
        scen_tbl_plant.lat, ...
        scen_tbl_plant.lon);                            
    params.market = params.D<=d_max;                  % market = counties within d_max miles of each plant             
    params.sigma = .01 + (.025/100)*min(params.D, ... % Shrink matrix: shrinkage (% of weight) to travel distances in D
        100*ones(size(params.D)));     
    params.T = t*params.D./(1-params.sigma);          % Transportation cost matrix: total transport cost for distances in D
    
    % NxN matrix of reporting regions: row n gives the counties that
    % county n uses to compute the regional average price used for AMAs
    params.reporting_regions = reportingRegions(tbl_county.region);

    % Define function handles to compute shares and marginal costs
    % These are small functions that need to be called from different
    % parent functions (initialPrice, computeNashEqm, etc.), so we define
    % these functions here and pass them as parameters
    params.shares = @(p) computeShares(p,params);
    params.marginal_cost = @(Q,K) params.c0 + (Q>K).*params.gamma.*(Q./K-1).^params.lambda;

    %% Import raw algorithm output to create results tables
    s_init = readmatrix(strcat('../Replication Package/raw_output/',name,'/s_init.csv')); % Initial shares used to compute contracting volumes for each county
    p_eqm = readmatrix(strcat('../Replication Package/raw_output/',name,'/p_eqm.csv'));   % Equilibrium prices, units of $/lb beef
    s_eqm = readmatrix(strcat('../Replication Package/raw_output/',name,'/s_eqm.csv'));   % Equilibrium shares
    mc_eqm = readmatrix(strcat('../Replication Package/raw_output/',name,'/mc_eqm.csv')); % Equilibrium firm-level marginal costs
    p_noc = readmatrix(strcat('../Replication Package/raw_output/',name,'/p_nc.csv'));    % No contract prices (for markdown decomp), units of $/lb beef
    p_ind = readmatrix(strcat('../Replication Package/raw_output/',name,'/p_ind.csv'));  % Independent plant prices (for markdown decomp), units of $/lb beef

    % Compute X: Contracted quantity in each county going to each plant.
    % X(j,n) = contracted quantity going from county n to plant j
    % Based initial shares and regional contract shares observed in data.
    X = (scen_tbl_county.contract_share)'.*s_init.*params.Q_n;

    % Marginal breakeven prices = highest price plants could offer without
    % losing money on the marginal head of cattle
    p_mbe = params.b-mc_eqm-params.T;

    % Unconstrained breakeven prices = highest price plants could offer if
    % there were no cap constraints (i.e. mc stayed at c0)
    p_be = params.b-params.c0-params.T;

    %% Save results

    % For mill price eqm - convert computed prices to prices net of
    % transportation costs for comparability with baseline eqm
    % First, save mill prices to results
    if scenario==8
        results.(name).p_eqm_mill = p_eqm*100*pct_beef;
        results.(name).p_noc_mill = p_noc*100*pct_beef;
        results.(name).p_ind_mill = p_ind*100*pct_beef;
        results.(name).p_mbe_mill = (params.b - mc_eqm)*100*pct_beef;
        results.(name).p_be_mill  = (params.b - params.c0)*100*pct_beef;

        p_eqm = p_eqm - params.T;
        p_noc = p_noc - params.T;
        p_ind = p_ind - params.T;
    end

    % AVERAGE PRICES AND MARKDOWN DECOMPOSITION
    % Markdown decomposition
    % Convert to $/cwt
    md_full    = (p_mbe - p_eqm) *100*pct_beef;     % Full markdown
    md_amas    = (p_noc - p_eqm) *100*pct_beef;     % Contracting markdown
    md_multi   = (p_ind - p_noc) *100*pct_beef;     % Multi-plant markdown
    md_spatial = (p_mbe - p_ind) *100*pct_beef;     % Spatial markdown
    cap_effect = (p_be - p_mbe)  *100*pct_beef;     % Capacity effect

    % County-level average prices and markdowns
    p_eqm_n      = sum(p_eqm.*s_eqm)*100*pct_beef;  % Average price (convert to $/cwt)
    md_full_n    = sum(md_full.*s_eqm);             % Full markdown
    md_amas_n    = sum(md_amas.*s_eqm);             % Contracting markdown
    md_multi_n   = sum(md_multi.*s_eqm);            % Multi-plant markdown
    md_spatial_n = sum(md_spatial.*s_eqm);          % Spatial markdown
    cap_effect_n = sum(cap_effect.*s_eqm);          % Capacity effect

    % National average prices and markdowns
    M_n = params.Q_n - sum(X); % Compute county-level spot market quantities for weighted avgs
    p_eqm_avg      = wavg(p_eqm_n,M_n);         % Average price
    md_full_avg    = wavg(md_full_n,M_n);       % Full markdown
    md_amas_avg    = wavg(md_amas_n,M_n);       % Contracting markdown
    md_multi_avg   = wavg(md_multi_n,M_n);      % Multi-plant markdown
    md_spatial_avg = wavg(md_spatial_n,M_n);    % Spatial markdown
    cap_effect_avg = wavg(cap_effect_n,M_n);    % Capacity effect
    
    % SAVE RESULTS
    % Save national averages for all scenarios
    results.(name).p_eqm_avg      = p_eqm_avg;
    results.(name).md_full_avg    = md_full_avg;
    results.(name).md_amas_avg    = md_amas_avg;
    results.(name).md_multi_avg   = md_multi_avg;
    results.(name).md_spatial_avg = md_spatial_avg;
    results.(name).cap_effect_avg = cap_effect_avg;

    % Save average distance cattle travel from feedlot to plant
    d_n = sum(params.D.*(s_eqm.*M_n + X))./params.Q_n;
    results.(name).d_avg = wavg(d_n,params.Q_n);

    % Save detailed raw results (this should be enough to compute any
    % model output of interest that isn't saved explicitly)
    results.(name).p_eqm      = p_eqm*100*pct_beef; % Convert to $/cwt (markdowns below are already in $/cwt)
    results.(name).md_full    = md_full;
    results.(name).md_spatial = md_spatial;
    results.(name).md_multi   = md_multi;
    results.(name).md_amas    = md_amas;
    results.(name).cap_effect = cap_effect;

    results.(name).s_eqm  = s_eqm;
    results.(name).mc_eqm = mc_eqm;
    results.(name).X      = X;
    results.(name).M_n    = M_n;

    % Save scenario params and environment
    results.(name).tbl_plant  = scen_tbl_plant;
    results.(name).params = params;

    %% Report results: Save tables and data for figures for manuscript
    
    % BASELINE and MILL PRICE EQUILIBRIA:
    % Save county-level results to table
    if scenario==1 || scenario==8
        % Generate table of county-level results
        state_name = tbl_county.state_name;
        fips = tbl_county.fips;
        FedCattle = tbl_county.fed_cattle;
        Price = p_eqm_n';
        Markdown = md_full_n';
        MarkdownSpatial = md_spatial_n';
        MarkdownMulti = md_multi_n';
        MarkdownContract = md_amas_n';
        CapEffect = cap_effect_n';

        % Store county-level results in table
        T = table(state_name,fips,FedCattle,Price,Markdown,MarkdownSpatial,MarkdownMulti,MarkdownContract,CapEffect);
        
        % For baseline eqm only, save table to results struct and to CSV
        if scenario==1
            results.(name).county_results = T;
            % Save table as CSV file
            % Used to make county-level maps - FIGURE 4 in manuscript
            writetable(T,"data/county_results.csv") % Save to data folder since this is used for ultimate output
        end

        % Summarize results for full U.S. and top 5 cattle producing states
        % TABLE 3 (Baseline) and TABLE 5 (Mill price) in manuscript
        results.(name).T_summary = eqmSummary(T,M_n);
    end
    
    % BASELINE EQUILIBRIUM and COLLUSION EQUILIBRIUM
    % Save firm-level results to table
    if scenario==1 || scenario==9
        % TABLE 4 in manuscript (Baseline) and results for supplementary
        % appendix on stability of collusion (Collusion eqm)
        if scenario==9
            results.(name).tbl_plant.firm_id = results.(name).tbl_plant.firm_id_save; % Change firm_id back to original values for proper firm aggregationg
        end
        results.(name).T_firm = firmSummary(p_eqm,s_eqm,results.(name),params);
    end

    % BASELINE AND COUNTERFACTUAL POLICY SCENARIOS
    % Save national avg results to table - TABLE 6 in manuscript
    if scenario<=7
        results.T_scenarios = scenarioResults(results,name);
    end

    % COLLUSION STABILITY
    % Report firm-level results when one firm unilaterally deviates from
    % collusive equilibrium
    % For supplemental appendix analyzing stability of collusion
    if scenario==9
        % Initialize table
        T = table();

        % Top 4 firms
        firms = {"JBS","Tyson","Cargill","National Beef"};
        
        % Loop over top 4 firms and compute firm-level results when they
        % deviate from collusion
        for i=1:4
            % Initialize price matrix to p_eqm
            p_dev = p_eqm;
    
            % Firm j deviates by increasing price by $0.10/cwt
            dev = 0.10/(100*pct_beef); % Convert to $/lb beef
            plants = scen_tbl_plant.firm_id_save==i;
            p_dev(plants,:) = p_dev(plants,:) + dev;
    
            % Compute shares after deviation
            s_dev = params.shares(p_dev);
            
            new_rows = firmSummary(p_dev,s_dev,results.(name),params);
            new_rows{:,"Deviator"} = firms{i};
            T = [T;new_rows];
        end
        results.(name).T_deviation = T;
    end
    


end

%% Functions
% Functions to create tables for reporting results

function T = eqmSummary(county_table,M_n)
% GENERATE TABLE 3 (Baseline)  AND TABLE 5 (Mill price) IN MANUSCRIPT
% Initialize table
T = table();

% Make M_n column vector to match table columns
M_n = M_n';

% Areas to report
areas = {"Full U.S.","Kansas","Nebraska","Texas","Colorado","Iowa"};

for i=1:length(areas)
    % Get counties to average over for the area
    if i==1
        counties = county_table.state_name~="";
    else
        counties = county_table.state_name==upper(areas{i});
    end
    
    % Get weighted average values for the area
    Area = areas{i};
    Price            = wavg(county_table.Price(counties),           M_n(counties));
    Markdown         = wavg(county_table.Markdown(counties),        M_n(counties));
    MarkdownSpatial  = wavg(county_table.MarkdownSpatial(counties), M_n(counties));
    MarkdownMulti    = wavg(county_table.MarkdownMulti(counties),   M_n(counties));
    MarkdownContract = wavg(county_table.MarkdownContract(counties),M_n(counties));
    CapEffect        = wavg(county_table.CapEffect(counties),       M_n(counties));
    MarkdownPct      = Markdown/Price;

    % Add to table
    new_row = table(Area,Price,Markdown,MarkdownPct,MarkdownSpatial,MarkdownMulti,MarkdownContract,CapEffect);
    T = [T;new_row];
end
end

function T = firmSummary(p,s,scen_results,params)
% GENERATE TABLE SUMMARIZING FIRM-LEVEL RESULTS
% Inputs: p and s: JxN matrices of plant- and county-specific prices and
% shares. These inputs are separate from the scen_results (where eqm p and
% s are stored) to allow function to compute firm-level results when one
% firm deviates from collusive equilibrium

% Initialize table
T = table();

% Get M_n and X as local variables
M_n = scen_results.M_n; % County-level spot market volume
X = scen_results.X;     % JxN matrix of contracted quantities

% Compute plant-level quantities
M_plant = sum(s.*M_n,2);        % Spot-market quantity
Q_plant = M_plant + sum(X,2);   % Total plant quantity

% Compute markdowns
% NOTE: Markdowns are saved in scen_results, but this function also needs
% to compute firm-level results when one firm deviates from collusion. So
% we need to compute markdowns for given p and s
mc = params.marginal_cost(Q_plant,params.K);
markdown = (params.b-mc-params.T) - p;

% Compute plant-level averages
p_plant = sum(p.*(s.*M_n),2)./M_plant*100*params.pct_beef;          % Plant-level average price (convert to $/cwt)
md_plant = sum(markdown.*(s.*M_n),2)./M_plant*100*params.pct_beef;  % Plant-level average markdown (convert to $/cwt)
d_plant = sum(params.D.*(s.*M_n + X),2)./Q_plant;                   % Plant-level average procurement distance
% NOTE: Markdown is defined in the spot-market, so only spot market
% quantites are used to compute weighted average. Average distance involves
% all cattle, so we compute weighted average over total plant quantity

% Compute plant-level profit
profit_plant = plantProfit(p,s,scen_results,params);

% "Firms" to report
firms = {"JBS","Tyson","Cargill","National Beef","Other","Industry"};

for i=1:length(firms)
    Firm = firms{i};

    % Get set of plants to average over
    if i<=4
        plants = scen_results.tbl_plant.firm_id==i;
    elseif i==5
        plants = scen_results.tbl_plant.firm_id > 4;
    else
        plants = scen_results.tbl_plant.firm_id > 0;
    end

    % Firm-level (weighted) averages:
    % Price, Markdown, Cap Utilization, Procurement Distance
    Price = wavg(p_plant(plants),M_plant(plants));
    Markdown = wavg(md_plant(plants),M_plant(plants));
    CapacityUtilization = sum(Q_plant(plants))/sum(params.K(plants));
    ProcurementDistance = wavg(d_plant(plants),Q_plant(plants));

    % Firm-level marginal costa  ratio to base mc (c0)
    AverageMarginalCost = mean(mc(plants))/params.c0;   % Simple average
    MaxMarginalCost = max(mc(plants))/params.c0;        % Max

    % Firm profits
    Profit = sum(profit_plant(plants));
    
    new_row = table(Firm,Price,Markdown,CapacityUtilization,ProcurementDistance,AverageMarginalCost,MaxMarginalCost,Profit);
    T = [T;new_row];
end
end

function T = scenarioResults(results,name)
% GENERATE TABLE WITH NATIONAL SUMMARY OF RESULTS FOR SCENARIOS --
% TABLE 6 IN MANUSCRIPT
% One row each for national average outcomes of baseline scenario and 
% 6 counterfactual scenarios
T = results.T_scenarios;
Scenario         = name;
Price            = results.(name).p_eqm_avg;
Markdown         = results.(name).md_full_avg;
MarkdownSpatial  = results.(name).md_spatial_avg;
MarkdownMulti    = results.(name).md_multi_avg;
MarkdownContract = results.(name).md_amas_avg;
CapEffect        = results.(name).cap_effect_avg;
MarkdownPct      = Markdown/Price;

new_row = table(Scenario,Price,MarkdownPct,Markdown,MarkdownSpatial,MarkdownMulti,MarkdownContract,CapEffect);
T = [T;new_row];
end

%% Save table files
writetable(results.base.T_summary,'results/Table_3.csv')
writetable(results.base.T_firm,'results/Table_4.csv')
writetable(results.mill.T_summary,'results/Table_5.csv')
writetable(results.T_scenarios,'results/Table_6.csv')
