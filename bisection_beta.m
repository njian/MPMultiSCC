function [f_n, x_ast_ast, SL_n, sd_n, beta_U, n_sim] = bisection_beta(n, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
epsilon = 0.1; % stopping precision for beta
beta_L = 0; % initial bounds for beta
beta_U = 1e4;
n_sim = 0;

count_beta = 0;
f_n = 1e4;

nShifts = size(shifts, 1);
x_trial = zeros(nAgentGroups, nShifts);      
for k = 1:n
    row = randi(nAgentGroups);
    col = randi(nShifts);
    x_trial(row, col) = x_trial(row, col) + 1;
end
x_ast = x_trial;
x_ast_ast = x_trial;

while beta_U - beta_L > epsilon
    z = max(2, log2(beta_U - beta_L));
    beta = beta_L + (beta_U - beta_L)/z; % new center point
    count_beta = count_beta + 1;
    fprintf('======================================================================================== \n');
    fprintf('Bisection iteration (%d) of beta, with n fixed at %d. \n', count_beta, n);
    fprintf('beta_L = %.2f, beta_U = %.2f. \n', beta_L, beta_U);

    if count_beta == 1 % initialize
        beta = beta_U;
        fprintf('Initialize: Try beta = beta_U first. \n');      
    else
       fprintf('beta = %.2f. \n', beta); 
    end
    
    % Local search for x given beta and n
    [f_beta, x_ast, SL_beta, sd_beta, reps] = localSearch_x(n, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
    %       [f_beta, x_ast, SL_beta, sd_beta, reps] = localSearch_varyN_x(n, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts); 
    n_sim = n_sim + reps;
    
    fprintf('======================================================================================== \n');
    if f_beta < f_n
        fprintf('beta step: IMPROVE f from %.2f to %.2f +/- %.2f. \n', f_n, f_beta, sd_beta);
        % update the best solution for any beta
        x_ast_ast = x_ast;
        fprintf('Number of agents: %d. \n', sum(sum(x_ast_ast)));
        f_n = f_beta;
        SL_n = SL_beta;
        sd_n = sd_beta;
    else
        x_ast = x_ast_ast;
        fprintf('beta step: failed: this solution %.2f >= last best %.2f. \n', f_beta, f_n);
    end

    % Bisection beta step: find beta s.t. SL = serviceLevelMin
    if count_beta == 1 && SL_beta < serviceLevelMin
        fprintf('Bisection of beta failed. SL_beta = %.2f, so n is not big enough for beta_U. \n', SL_beta);
        break
    end

    if SL_beta < serviceLevelMin
        beta_L = beta;
    else
        beta_U = beta;
    end

end % end bisection search for beta

end

