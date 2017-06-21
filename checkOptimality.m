% Check optimality
clear;
load('C:\Users\nj227\Dropbox\Research\CallCenter\MPMultiSCC\outputs\special\util\large_seed1.mat');
[nShifts, ~] = size(shifts);
[f, SL, ~, ~, ~, ~, avgServerUtil0] = MultiSkillPickedCalls(x_opt, beta, 30, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts);
avgServerUtilByGroup0 = nanmean(avgServerUtil0,2);
improvedComponents = [0, 0 f, SL, avgServerUtilByGroup0(1), avgServerUtilByGroup0(2)];

for g = 1:nAgentGroups
    for s = 1:nShifts
        if x_opt(g,s) >= 1
            x = x_opt;
            x(g,s) = x(g,s) - 1;
            [f1, SL1, ~, ~, ~, ~, avgServerUtil1] = MultiSkillPickedCalls(x, beta, 30, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts);
            if SL1 > 0.8
                avgServerUtilByGroup1 = nanmean(avgServerUtil1,2);
                improvedComponents = [improvedComponents; [g, s, f1, SL1, avgServerUtilByGroup1(1), avgServerUtilByGroup1(2)]];
            end
        end 
    end
end
save('C:\Users\nj227\Dropbox\Research\CallCenter\MPMultiSCC\outputs\special\util\large_seed1_optimality.mat');