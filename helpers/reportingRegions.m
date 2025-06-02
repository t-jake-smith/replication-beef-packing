function reporting_regions = reportingRegions(region)
% reporting_regions: For each county (n) in the vector of the region 
% each of the N counties belongs to (region), find 1xN vector indicating
% which counties are in the area that county n uses for its average price.

% The areas that each region uses for average price are:
% REGION 1 (TX/OK/NM): TX/OK/NM
% REGION 2 (KS): KS
% REGION 3 (NE): NE
% REGION 4 (CO): 5 area (regions 1:5)
% REGION 5 (IA/MN/MO): 5 area (regions 1:5)
% REGIONS 6-10 (all other regions): National average price

% Number of counties
N = length(region);
% Make region a row vector
region = region';
% NxN reporting_regions matrix. Row n gives the counties used to
% compute p_bar (average price) for county n
reporting_regions = zeros(N);

for n=1:N
    county_region = region(n);
    
    % TX/OK/NM, KS, and NE regions: compute p_bar for their own specific region
    if county_region==1 || county_region==2 || county_region==3
        reporting_regions(n,:) = region==county_region;
    
    % IA/MN/MO and CO: compute p_bar for the 5 Area greater region
    elseif county_region==4 || county_region==5
        reporting_regions(n,:) = ismember(region,1:5);
    
    % All other regions: compute p_bar at national level
    else
        reporting_regions(n,:) = region>0;
    end
end
reporting_regions = logical(reporting_regions);
end