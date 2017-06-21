clear;
path = 'C:\Users\nj227\Dropbox\Research\CallCenter\MPMultiSCC\outputs\special\search_n\combinationSearch\small\';
filename = '05_24_00_44_51';
load([path, filename, '.mat']);

[nGroup, nShift] = size(x_opt);
obj_m = zeros(nGroup, nShift);
SL_m = zeros(nGroup, nShift);
obj_p = zeros(nGroup, nShift);
SL_p = zeros(nGroup, nShift);
potential_add = 0;
potential_remove = 0;
runlength = 10;

tic;
for group = 1:nGroup
    for shift = 1:nShift
        x_trial = x_opt;
        if x_trial(group, shift) > 0
            x_trial(group, shift) = x_trial(group, shift) - 1;
            [obj_m(group, shift), SL_m(group, shift)] = MultiSkillPickedCalls(x_trial, beta, ceil(runlength), seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
        else
            obj_m(group, shift) = NaN;
            SL_m(group, shift) = NaN;
        end
        
        if SL_m(group, shift) >= serviceLevelMin && obj_m(group, shift) > f_opt
            potential_remove = potential_remove + 1;
        end
        
        x_trial = x_opt;
        x_trial(group, shift) = x_trial(group, shift) + 1;
        [obj_p(group, shift), SL_p(group, shift)] = MultiSkillPickedCalls(x_trial, beta, ceil(runlength), seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
    
        if SL_p(group, shift) >= serviceLevelMin && obj_p(group, shift) > f_opt
            potential_add = potential_add + 1; 
        end
    end
end
grad_m = obj_m - f_opt;
grad_p = obj_p - f_opt;
fprintf('Number of potential adds: %d. \n', potential_add);
fprintf('Number of potential removals: %d. \n', potential_remove);
t = toc

save([path, filename, 'CheckOptimality']);
