function [f_beta, x_ast, SL_beta, sd_beta, n_sim] = localSearch_x(n, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts)
% Given n and beta, local search for x.
nShifts = size(shifts, 1);

% Search Parameters
nChangeMax = 32; % min(32, 2^floor(log2(nAgentGroups*nShifts/10))); % maximum number of swaps in x when generating trial solution
r = 5; % local random search among this many largest gradient components
max_fail = 5; % maximum number of consecutive fails in local search for x

count_x = 0;
nChange = nChangeMax; 
count_failed_x = 0;
n_sim = 0;
% x_trial = x_trial_in;

% Initialize
x_trial = zeros(nAgentGroups, nShifts);    
for k = 1:n
    row = randi(nAgentGroups);
    col = randi(nShifts);
    x_trial(row, col) = x_trial(row, col) + 1;
end
[f_beta_x, SL_beta_x, sd_x, forwardGradient_x, backwardGradient_x] = MultiSkillPickedCalls(x_trial, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts);
f_beta = f_beta_x;
SL_beta = SL_beta_x;
sd_beta = sd_x;
forwardGradient = forwardGradient_x;
backwardGradient = backwardGradient_x;
backwardGradient(x_trial==0) = -Inf;

% f_beta = f_beta_x;
x_ast = x_trial;
% SL_beta = SL_beta_x;
% sd_beta = sd_beta_x;
STOP = 0;
% stop when either consecutive fails exceed tolerance, or there's
% no room for improvement
while count_failed_x < max_fail && STOP == 0
    count_x = count_x + 1;
    fprintf('--------------------------------------------------------------------- \n');
    fprintf('Local search iteration (%d) of x. \n', count_x);
    fprintf('n = %d, beta = %.2f, nChange = %d. \n', sum(sum(x_trial)), beta, nChange);

    % Generate trial solution x: ramndomly choose among the largest
    % components in forward and backward gradients.
    % Each group*shift will only be changed at most once.
    % WIP: if more, need to deal with below 0 in mMax case.      
    % backwardGradient = - forwardGradient;
    if count_failed_x == 0
        [~, IndForwardSave] = sort(forwardGradient(:),'descend');
        [~, IndBackwardSave] = sort(backwardGradient(:),'descend');
    end
    IndForward = IndForwardSave;
    IndBackward = IndBackwardSave;
    % the positive components in forwardGradient and non-negative
    % components in backwardGradient indicate room for improvement
    npForward = sum(forwardGradient(:)>=0);
    npBackward = sum(backwardGradient(:)>=0);
    if npForward <=0 && npBackward <= 0
        STOP = 3; % local maximum
    elseif npForward <= 0 
        STOP = 1; % all coordinates want to decrease
    elseif npBackward <= 0
        STOP = 2; % all coordinates want to increase
    end
    fprintf('Number of positive forward and backward gradients: %d and %d \n', npForward, npBackward);

    % From Cezik 2008: adaptive step size to avoid non-concave regions
    if SL_beta < 0.5
        d = 3;
    elseif SL_beta < 0.65
        d = 2;
    else
        d = 1;
    end
    
    changeTolerancePlus = d * ones(nAgentGroups, nShifts);
    changeToleranceMinus = d * ones(nAgentGroups, nShifts);
    countChange = 0;
    while countChange < nChange && npForward > 0 && npBackward > 0 && STOP == 0 
        % search radius = min(r, number of positive components)
        rF = min(r, npForward);
        rB = min(r, npBackward);
        
        % Add agents
        randIndexF = randi(rF);
        temp = IndForward(randIndexF);
        [pRow, pCol] = ind2sub([nAgentGroups, nShifts], temp);
        x_trial(pRow, pCol) = x_trial(pRow, pCol) + 1;
        changeTolerancePlus(pRow, pCol) = changeTolerancePlus(pRow, pCol) - 1;
        if changeTolerancePlus(pRow, pCol) <= 0
            IndForward(randIndexF) = []; % remove the chosen index
            npForward = npForward - 1;
        end

        % Remove agents
        randIndexB = randi(rB);
        temp = IndBackward(randIndexB);
        [mRow, mCol] = ind2sub([nAgentGroups, nShifts], temp);

        x_trial(mRow, mCol) = x_trial(mRow, mCol) - 1; 
        changeToleranceMinus(mRow, mCol) = changeToleranceMinus(mRow, mCol) - 1;
        if (changeToleranceMinus(mRow, mCol) <= 0 || x_trial(mRow, mCol) <= 0)
            IndBackward(randIndexB) = [];
            npBackward = npBackward - 1;
        end
        countChange = countChange + 1;        
    end
    
    if countChange == 0
       fprintf('x step: no trial solution can be generated. \n'); 
       STOP = 4;
    end
    
    if STOP == 0 && countChange > 0
        [f_beta_x, SL_beta_x, sd_x, forwardGradient_x, backwardGradient_x] = MultiSkillPickedCalls(x_trial, beta, ceil(runlength), seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts);
        fprintf('Trial solution obj = %.2f, SL = %.2f. \n', f_beta_x, SL_beta_x);        
        n_sim = n_sim + runlength;

        if f_beta_x > f_beta
            fprintf('x step: IMPROVE f_beta from %.2f to %.2f +/- %.2f. \n', f_beta, f_beta_x, sd_x);
%             fprintf('x step: n = %d. \n', sum(sum(x_trial)));
            count_failed_x = 0;

            % update last best solution for this beta
            x_ast = x_trial;
            f_beta = f_beta_x;
            sd_beta = sd_x;
            SL_beta = SL_beta_x;
            forwardGradient = forwardGradient_x;
            backwardGradient = backwardGradient_x;
            backwardGradient(x_ast==0) = -Inf; % remove impossible reductions in x
        else
            x_trial = x_ast;
            count_failed_x = count_failed_x + 1;
            fprintf('x step: failed for %d times: this solution %.2f < last best %.2f. \n', count_failed_x, f_beta_x, f_beta);
            % update stepsize
            if count_failed_x == max_fail && nChange > 1
                nChange = ceil(nChange/2);
                count_failed_x = 0;
            end

        end
    elseif STOP ~= 4
        % First parse found no room for improvement. Need to
        % increase/decrease n.
        fprintf('x step: no room for improvement. ');
        if STOP == 3
            fprintf('Seems like a local optimal. \n');
        elseif STOP == 2
            fprintf('n is too small for this level of beta. \n');
        elseif STOP == 1
            fprintf('n is too big for this level of beta. \n');
        end
    end
end % end local search for x

end

