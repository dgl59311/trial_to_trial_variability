function [g_avg] = delete_bl(g_avg)
trial_starts = find(g_avg.times == 0);
g_avg.times = g_avg.times(trial_starts:end);
g_avg.data = g_avg.data(:, trial_starts:end);
g_avg.pnts = length(g_avg.times);
g_avg.xmin = g_avg.times(1);
end
