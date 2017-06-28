function [ x_opt, beta, history_n, history_all, n_sim, n_beta] = goldenSection_n( runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime, costByGroup)
% Use golden section to find n with minimal f_n = max_n f_n_beta.
golden_ratio = 2/(1 + sqrt(5));

% EXCLUSIVE BOUNDS: evaluating all n from n_L + 1 to n_U - 1.
if isempty(shifts)
    shiftsCoverage = 1;
else
    shiftsCoverage = ceil( max(shifts(:,5))/min(shifts(:,6)) );
end

% LB: treating as single-skill and use fastest agent group
n_L = solve_ErlangC( max(sum(arrivalRates,1))/size(arrivalRates,2), min(meanST(:)), SLtime ); %mean(arrivalRates(i,:)) * min(meanST(i,:));
n_U = 0;
for i = 1:nCallTypes  
   % UB: each callType separate call center, using slowest agent group in
   % each
   c2 = solve_ErlangC( max(arrivalRates(i,:)), max(meanST(i,:)), SLtime );
   n_U = n_U + c2; 
end
n_L = floor(n_L*shiftsCoverage); 
n_U = ceil(n_U*shiftsCoverage);
% Test if the bounds are appropriate
[nShifts, ~] = size(shifts);
% stepSize = max(nAgentGroups * nShifts * 0.1, 5);
% x_trial = evenly_spread( n_U, nAgentGroups, nShifts );
% [~, ~, ~, forwardGradient, ~] = MultiSkillPickedCallsPots(x_trial, 1e3, 1, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime);
% while sum(sum(forwardGradient>0)) > 0.8 * nAgentGroups * nShifts
%     n_U = n_U + stepSize;
%     x_trial = evenly_spread( n_U, nAgentGroups, nShifts );
%     [~, ~, ~, forwardGradient, ~] = MultiSkillPickedCallsPots(x_trial, 1e3, 1, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime);
% end
% x_trial = evenly_spread( n_L, nAgentGroups, nShifts );
% [~, ~, ~, ~, backwardGradient] = MultiSkillPickedCallsPots(x_trial, 1e3, 1, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime);
% while n_L >= stepSize && sum(sum(backwardGradient>0)) > 0.8 * nAgentGroups * nShifts
%     n_L = n_L - stepSize;
%     x_trial = evenly_spread( n_L, nAgentGroups, nShifts );
%     [~, ~, ~, ~, backwardGradient] = MultiSkillPickedCallsPots(x_trial, 1e3, 1, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime);
% end

% [~, ~, SL] = localSearch_x_aggressive(n_U, 1e3, 1, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime, costByGroup);
% while SL < 0.8
%     n_U = n_U + stepSize;
%     [~, ~, SL] = localSearch_x_aggressive(n_U, 1e3, 1, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime, costByGroup);
% end
% [~, ~, SL] = localSearch_x_aggressive(n_L, 1, 1, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime, costByGroup);
% while SL > 0.8 && n_L > stepSize
%     n_L = n_L - stepSize;
%     [~, ~, SL] = localSearch_x_aggressive(n_L, 1, 1, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime, costByGroup);
% end

if n_U - n_L <= 1
    fprintf('Error: Invalid n_L and n_U bounds. Bounds are exclusive: if want n = 1, one should set n_L = 0, n_U = 3. \n');
end

history_n = [];
history_all = [];
f = -1e4; % best solution so far
count_n = 0;
n_sim = 0;
n_beta = 0;
x_opt = NaN;
beta = NaN;

reuse_rho = 0;
reuse_lambda = 0;
while n_U - n_L > 1
    count_n = count_n + 1;
    fprintf('\n\n\n*************************************************************************************************************************** \n');
    fprintf('*************************************************************************************************************************** \n');
    fprintf('Golden section iteration (%d) of n: n_L = %i, n_U = %i. \n', count_n, n_L, n_U);
    fprintf('n_sim = %d. \n', n_sim);
    
    aggressiveSearch = 1;
%     if n_U - n_L <= 10 % approximately 4 iterations to go
%         aggressiveSearch = 1;
%     end  
    fprintf('Using aggressiveSearch? %d \n', aggressiveSearch);
    
    % First try to reuse n_rho or n_lambda points utilizing golden ratio.
    if reuse_lambda == 1
        n_rho = n_lambda;
        x_rho = x_lambda;
        f_rho = f_lambda;
        beta_rho = beta_lambda;
        SL_rho = SL_lambda;
        sd_rho = sd_lambda;
    end
    if reuse_rho == 1
        n_lambda = n_rho;
        x_lambda = x_rho;
        f_lambda = f_rho;
        beta_lambda = beta_rho;
        SL_lambda = SL_rho;
        sd_lambda = sd_rho;
    end
    % Then evaluate the missing n_rho or n_lambda.
    skip_rest = 0;
    if reuse_lambda == 0 % cannot reuse old lambda for new rho
        if count_n == 1 || n_rho ~= round((1 - golden_ratio) * n_L + golden_ratio * n_U) % if not new n_rho = new n_lambda = old n_rho
            n_rho = round((1 - golden_ratio) * n_L + golden_ratio * n_U);
            fprintf('*************************************************************************************************************************** \n');
            fprintf('Evaluating: n_rho = %i. \n', n_rho);
            [f_rho, x_rho, SL_rho, sd_rho, beta_rho, reps, fAll_beta, count_beta] = bisection_beta(n_rho, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime, aggressiveSearch, costByGroup);
            history_all = [history_all, fAll_beta];
            n_sim = n_sim + reps;
            n_beta = n_beta + count_beta;
        end            
        % If SL is too low, increase n directly to the larger point n_rho, 
        % and skip evaluating the smaller point n_lambda. 
        if SL_rho < serviceLevelMin
            skip_rest = 1;
        end
        fprintf('*************************************************************************************************************************** \n');
        fprintf('Concluding: n_rho = %i, f_rho = %.2f, SL_rho = %.2f. skip_rest = %d. \n', n_rho, f_rho, SL_rho, skip_rest);
    end
    
    if skip_rest == 0
        if reuse_rho == 0 % cannot reuse old rho for new lambda
            n_lambda = round(golden_ratio * n_L + (1 - golden_ratio) * n_U);
            if n_lambda ~= n_rho
                fprintf('*************************************************************************************************************************** \n');
                fprintf('Evaluating: n_lambda = %d. \n', n_lambda);
                [f_lambda, x_lambda, SL_lambda, sd_lambda, beta_lambda, reps, fAll_beta, count_beta] = bisection_beta(n_lambda, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime, aggressiveSearch, costByGroup);
                history_all = [history_all, fAll_beta];
                n_sim = n_sim + reps;
                n_beta = n_beta + count_beta;
            else
                n_lambda = n_rho;
                x_lambda = x_rho;
                f_lambda = f_rho;
                beta_lambda = beta_rho;
                SL_lambda = SL_rho;
                sd_lambda = sd_rho;
            end
        end
        fprintf('*************************************************************************************************************************** \n');
        fprintf('Concluding: n_lambda = %i, f_lambda = %.2f, SL_lambda = %.2f. \n', n_lambda, f_lambda, SL_lambda);
    end
    fprintf('*************************************************************************************************************************** \n');
    fprintf('Finished golden section iteration (%d) of n: n_L = %i, n_U = %i. \n', count_n, n_L, n_U);
    
    % Increase n if SL is too low
    if SL_rho < serviceLevelMin % skip_rest = 1
        n_L = n_rho;
        f_n = f_rho;
        SL_n = SL_rho;
        sd_n = sd_rho;
        beta_n = beta_rho;
        x_n = x_rho;
        reuse_rho = 0; reuse_lambda = 0;
    % Increase n if SL is too low
    elseif SL_lambda < serviceLevelMin
        n_L = n_lambda;
        f_n = f_lambda;
        SL_n = SL_lambda;
        sd_n = sd_lambda;
        beta_n = beta_lambda;
        x_n = x_lambda;
        reuse_rho = 1; reuse_lambda = 0;
    elseif n_lambda == n_rho % when n_U - n_L = 4 or 2
        n_U = n_U - 1;
        f_n = f_lambda;
        SL_n = SL_lambda;
        sd_n = sd_lambda;
        beta_n = beta_lambda;
        x_n = x_lambda;
        reuse_lambda = 1; reuse_rho = 0;
    % If both SL satisfy SLMin, update bounds like golden ratio search
    elseif f_lambda < f_rho
        n_L = n_lambda;
        f_n = f_rho;
        SL_n = SL_rho;
        sd_n = sd_rho;
        beta_n = beta_rho;
        x_n = x_rho;
        reuse_rho = 1; reuse_lambda = 0;
    elseif f_lambda >= f_rho
        n_U = n_rho;
        f_n = f_lambda;
        SL_n = SL_lambda;
        sd_n = sd_lambda;
        beta_n = beta_lambda;
        x_n = x_lambda;
        reuse_lambda = 1; reuse_rho = 0;
    end        
    
    cost_n = beta * (SL_n - serviceLevelMin) - f_n;
    history_n = [history_n, [sum(sum(x_n)); n_sim; f_n; sd_n; cost_n; SL_n; beta_n]];
    
    % Record the best solution so far
    if f_n >= f        
        fprintf('n step: improve f from %.2f to %.2f +/- %.2f. SL = %.2f, beta = %.2f. \n', f, f_n, sd_n, SL_n, beta_n);
        f = f_n;
        x_opt = x_n;
        beta = beta_n;
    else
        fprintf('n step: failed. This solution %.2f < last best %.2f. \n', f_n, f);
    end
    
    if skip_rest == 0
        fprintf('Making n step: n_lambda = %i, f_lambda = %.2f, n_rho = %i, f_rho = %.2f. \n', n_lambda, f_lambda, n_rho, f_rho);
    else
        fprintf('n_rho = %i, f_rho = %.2f is not feasible. Increasing n.\n', n_rho, f_rho);
    end
    fprintf('*************************************************************************************************************************** \n');
    
end % end goldensection search for n

end