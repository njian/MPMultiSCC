function [ x_opt, beta, history_n, history_all, n_sim] = goldenSection_n( runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts)
% Use golden section to find n with minimal f_n = max_n f_n_beta.
golden_ratio = 2/(1 + sqrt(5));

% EXCLUSIVE BOUNDS: evaluating all n from n_L + 1 to n_U - 1.
n_L = 0; % initial bounds for #agents
% nCallTypes, nAgentGroups
[minArrivalRate, temp] = min(arrivalRates(:));
[minCallType, ~] = ind2sub(size(arrivalRates), temp);
n_U = ceil( minArrivalRate * nCallTypes * max(meanST(minCallType, :))  / 15); % required by 15-min max workload required
if n_U - n_L <= 1
    fprintf('Error: Invalid n_L and n_U bounds. Bounds are exclusive: if want n = 1, one should set n_L = 0, n_U = 3. \n');
end

history_n = [];
history_all = [];
f = -1e4; % best solution so far
count_n = 0;
n_sim = 0;
x_opt = NaN;
beta = NaN;

reuse_rho = 0;
reuse_lambda = 0;
while n_U - n_L > 1
    count_n = count_n + 1;
    fprintf('\n\n\n*************************************************************************************************************************** \n');
    fprintf('*************************************************************************************************************************** \n');
    fprintf('Golden section iteration (%d) of n: n_L = %d, n_U = %d. \n', count_n, n_L, n_U);
    fprintf('n_sim = %d. \n', n_sim);
    
    aggressiveSearch = 0;
    if n_U - n_L <= 10 % approximately 4 iterations to go
        aggressiveSearch = 1;
    end  
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
            fprintf('Evaluating: n_rho = %d. \n', n_rho);
            [f_rho, x_rho, SL_rho, sd_rho, beta_rho, reps, fAll_beta] = bisection_beta(n_rho, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, aggressiveSearch);
            history_all = [history_all, fAll_beta];
            n_sim = n_sim + reps;
        end            
        % If SL is too low, increase n directly to the larger point n_rho, 
        % and skip evaluating the smaller point n_lambda. 
        if SL_rho < serviceLevelMin
            skip_rest = 1;
        end
        fprintf('*************************************************************************************************************************** \n');
        fprintf('Concluding: n_rho = %d, f_rho = %.2f, SL_rho = %.2f. skip_rest = %d. \n', n_rho, f_rho, SL_rho, skip_rest);
    end
    
    if skip_rest == 0
        if reuse_rho == 0 % cannot reuse old rho for new lambda
            n_lambda = round(golden_ratio * n_L + (1 - golden_ratio) * n_U);
            if n_lambda ~= n_rho
                fprintf('*************************************************************************************************************************** \n');
                fprintf('Evaluating: n_lambda = %d. \n', n_lambda);
                [f_lambda, x_lambda, SL_lambda, sd_lambda, beta_lambda, reps, fAll_beta] = bisection_beta(n_lambda, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, aggressiveSearch);
                history_all = [history_all, fAll_beta];
                n_sim = n_sim + reps;
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
        fprintf('Concluding: n_lambda = %d, f_lambda = %.2f, SL_lambda = %.2f. \n', n_lambda, f_lambda, SL_lambda);
    end
    fprintf('*************************************************************************************************************************** \n');
    fprintf('Finished golden section iteration (%d) of n: n_L = %d, n_U = %d. \n', count_n, n_L, n_U);
    
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
        fprintf('Making n step: n_lambda = %d, f_lambda = %.2f, n_rho = %d, f_rho = %.2f. \n', n_lambda, f_lambda, n_rho, f_rho);
    else
        fprintf('n_rho = %d, f_rho = %.2f is not feasible. Increasing n.\n', n_rho, f_rho);
    end
    fprintf('*************************************************************************************************************************** \n');
    
end % end goldensection search for n

end