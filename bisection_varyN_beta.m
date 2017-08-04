function [f_n, x_ast_ast, SL_n, sd_n, beta_U, n_sim, fAll_beta] = bisection_varyN_beta(runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
epsilon = 0.1; % stopping precision for beta
beta_L = 0; % initial bounds for beta
beta_U = 1e3;
n_sim = 0;

count_beta = 0;
f_n = 1e4;
x_ast_ast = NaN;
SL_n = 0;
sd_n = 0;
fAll_beta = [];

% nShifts = size(shifts, 1);
% x_trial = zeros(nAgentGroups, nShifts);      
% for k = 1:n
%     row = randi(nAgentGroups);
%     col = randi(nShifts);
%     x_trial(row, col) = x_trial(row, col) + 1;
% end
% x_ast = x_trial;
% x_ast_ast = x_trial;

while beta_U - beta_L > epsilon
    z = max(2, log2(beta_U - beta_L));
    beta = beta_L + (beta_U - beta_L)/z; % new center point
    count_beta = count_beta + 1;
    fprintf('======================================================================================== \n');
    fprintf('Bisection iteration (%d) of beta. \n', count_beta);
    fprintf('beta_L = %.2f, beta_U = %.2f. \n', beta_L, beta_U);

    fprintf('beta = %.2f. \n', beta); 
    
    % Local search for x given beta and n
    [f_beta, x_ast, SL_beta, sd_beta, reps] = localSearch_varyN_x(x_ast_ast, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, SLtime, arrivalRates, meanST, R, Route, shifts); 
    n_sim = n_sim + reps;
    
    cost_beta = beta * (SL_beta - serviceLevelMin) - f_beta;
    fAll_beta = [fAll_beta, [sum(sum(x_ast)); n_sim; f_beta; sd_beta; cost_beta; SL_beta; beta]];
    fprintf('======================================================================================== \n');
    if f_beta < f_n %&& SL_beta >= serviceLevelMin
        fprintf('beta step: IMPROVE f from %.2f to %.2f +/- %.2f. \n', f_n, f_beta, sd_beta);
        % update the best solution for any beta
        x_ast_ast = x_ast;
        f_n = f_beta;
        SL_n = SL_beta;
        sd_n = sd_beta;
%         cost_n = beta * (SL_n - serviceLevelMin) - f_n;
%         fAll_beta = [fAll_beta, [sum(sum(x_ast_ast)); n_sim; f_n; sd_n; cost_n; SL_n; beta]];
        fprintf('Number of agents = %d. SL = %.2f. \n', sum(sum(x_ast_ast)), SL_n);
%     elseif SL_beta < serviceLevelMin
%         fprintf('beta step: Infeasible SL, increasing beta_L from %.2f to %.2f. \n', beta_L, beta);
    else
        fprintf('beta step: failed: this solution %.2f >= last best %.2f. Last SL = %.2f. \n', f_beta, f_n, SL_n);
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
    
    fprintf('======================================================================================== \n');
end % end bisection search for beta

end

