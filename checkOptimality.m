% Check optimality
clear;
load('C:\Users\nj227\Dropbox\Research\CallCenter\MPMultiSCC\outputs\special\util\Simopt\small_seed1.mat');
if isempty(shifts)
    nShifts = 1;
else
    [nShifts, ~] = size(shifts);
end
[f, SL, ~, ~, ~, ~, avgServerUtil0] = MultiSkillPickedCallsScheduling(x_opt, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime);
avgServerUtilByGroup0 = nanmean(avgServerUtil0,2);
improvedComponents_m = [0, 0 f, SL, avgServerUtilByGroup0(1), avgServerUtilByGroup0(2)];
improvedComponents_p = [0, 0 f, SL, avgServerUtilByGroup0(1), avgServerUtilByGroup0(2)];

for g = 1:nAgentGroups
    for s = 1:nShifts
        if x_opt(g,s) >= 1
            x1 = x_opt;
            x1(g,s) = x1(g,s) - 1;
            [f1, SL1, ~, ~, ~, ~, avgServerUtil1] = MultiSkillPickedCallsScheduling(x1, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime);
            if SL1 > 0.8 && f1 > f
                avgServerUtilByGroup1 = nanmean(avgServerUtil1,2);
                improvedComponents_m = [improvedComponents_m; [g, s, f1, SL1, avgServerUtilByGroup1(1), avgServerUtilByGroup1(2)]];
            end
        end 
    end
end

for g = 1:nAgentGroups
    for s = 1:nShifts
        x2 = x_opt;
        x2(g,s) = x2(g,s) + 1;
        [f2, SL2, ~, ~, ~, ~, avgServerUtil2] = MultiSkillPickedCallsScheduling(x2, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime);
        if SL2 > 0.8 && f2 > f
            avgServerUtilByGroup2 = nanmean(avgServerUtil2,2);
            improvedComponents_p = [improvedComponents_p; [g, s, f2, SL2, avgServerUtilByGroup2(1), avgServerUtilByGroup2(2)]];
        end
    end
end
% save('C:\Users\nj227\Dropbox\Research\CallCenter\MPMultiSCC\outputs\special\util\Pots_seed1_optimality.mat');