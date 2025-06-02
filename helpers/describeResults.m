function description = describeResults(array)
% describeResults: Get an array of prices, markdowns, or anything else, and
% return the min, max, average, median, 5th percentile, and 95th
% percentile, 1st percentile, and 99th percentile
% (These values are useful for deciding what range to use for mapping
% results)

description.min = min(array);
description.max = max(array);
description.median = median(array);
description.mean = mean(array);
description.p1 = prctile(array,1);
description.p5 = prctile(array,5);
description.p10 = prctile(array,10);
description.p90 = prctile(array,90);
description.p95 = prctile(array,95);
description.p99 = prctile(array,99);
end