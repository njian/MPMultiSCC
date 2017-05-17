function [f_beta, x_ast, SL_beta, sd_beta, n_sim] = localSearch_x(n, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts)
% Given n and beta, local search for x.
nShifts = size(shifts, 1);

% Search Parameters
nChangeMax = 2^floor(log2(nAgentGroups*nShifts/10)); % maximum number of swaps in x when generating trial solution
r = 5; % local random search among this many largest gradient components
max_fail = 10; % maximum number of consecutive fails in local search for x

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
[f_beta, SL_beta, sd_beta, forwardGradient, backwardGradient] = MultiSkillPickedCalls(x_trial, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
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
    [~, IndForward] = sort(forwardGradient(:),'descend');
    [~, IndBackward] = sort(backwardGradient(:),'descend');
    % the positive components in forwardGradient and non-negative
    % components in backwardGradient indicate room for improvement
    npForward = sum(forwardGradient(:)>=0);
    npBackward = sum(backwardGradient(:)>=0);
    if npForward <=0 && npBackward <= 0
        STOP = 1; % local maximum
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
%         nChange = min(nChange, d*(npForward+npBackward)); % e.g. if there are only 2 positive gradient components, make nChange <= 2
    changeTolerancePlus = d * ones(nAgentGroups, nShifts);
    changeToleranceMinus = d * ones(nAgentGroups, nShifts);
    countChange = 0;
    
    % search radius = min(r, number of positive components)
    [f1, f2] = ind2sub(size(forwardGradient), IndForward(1)); % index of largest forwardGradient
    rB = min(r, sum(sum(backwardGradient > -forwardGradient(f1, f2))));
    [f3, f4] = ind2sub(size(backwardGradient), IndBackward(1)); % index of largest forwardGradient
    rF = min(r, sum(sum(forwardGradient > -backwardGradient(f3, f4))));
    while countChange < nChange && rF > 0 && rB > 0 && STOP == 0 %(npForward > 0 || npBackward > 0)
        % Add agents
        randIndex = randi(rF);
        [pRow, pCol] = ind2sub(size(forwardGradient), IndForward(randIndex));
        % Remove agents
        randIndex = randi(rB);
        [mRow, mCol] = ind2sub(size(backwardGradient), IndBackward(randIndex));

        if forwardGradient(pRow, pCol) + backwardGradient(mRow, mCol) > 0
            x_trial(pRow, pCol) = x_trial(pRow, pCol) + 1;
            changeTolerancePlus(pRow, pCol) = changeTolerancePlus(pRow, pCol) - 1;
            if changeTolerancePlus(pRow, pCol) <= 0
                IndForward(randIndex) = []; % remove the chosen index
%                 npForward = npForward - 1;
            end
            x_trial(mRow, mCol) = x_trial(mRow, mCol) - 1; 
            changeToleranceMinus(mRow, mCol) = changeToleranceMinus(mRow, mCol) - 1;
            if (changeToleranceMinus(mRow, mCol) <= 0 || x_trial(mRow, mCol) <= 0)
                IndBackward(randIndex) = [];
%                 npBackward = npBackward - 1;
            end
            countChange = countChange + 1;        
        end
        
        % search radius = min(r, number of positive components)
        [f1, f2] = ind2sub(size(forwardGradient), IndForward(1)); % index of largest forwardGradient
        rB = min(r, sum(sum(backwardGradient > -forwardGradient(f1, f2))));
        [f3, f4] = ind2sub(size(backwardGradient), IndBackward(1)); % index of largest forwardGradient
        rF = min(r, sum(sum(forwardGradient > -backwardGradient(f3, f4))));
    end
    
    if STOP == 0 %~(STOP ~= 0 && count_x == 1)
        [f_beta_x, SL_beta_x, sd_x, forwardGradient_x, backwardGradient_x] = MultiSkillPickedCalls(x_trial, beta, ceil(runlength), seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
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
    else
        % First parse found no room for improvement. Need to
        % increase/decrease n.
        fprintf('x step: no room for improvement. Seems like a local optimal. \n');
    end
end % end local search for x

end

