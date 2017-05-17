function [f_beta, x_ast, SL_beta, sd_beta, n_sim] = localSearch_varyN_x(n, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts)
% Given n and beta, local search for x.
nShifts = size(shifts, 1);

% Search Parameters
nChangeMax = 1; %2^floor(log2(nAgentGroups*nShifts/10)); % maximum number of swaps in x when generating trial solution
r = 5; % local random search among this many largest gradient components
max_fail = 10; % maximum number of consecutive fails in local search for x

count_x = 0;
nChange = nChangeMax; 
count_failed_x = 0;
n_sim = 0;
% x_trial = x_trial_in;

% Initialize
% n_in_group = sum(x_trial_in, 2);
% x_trial = zeros(nAgentGroups, nShifts);
% for i = 1:nAgentGroups
%     while n_in_group(i) > 0
%         s = randi(nShifts);
%         x_trial(i,s) = x_trial(i,s) + 1;
%         n_in_group(i) = n_in_group(i) - 1;
%     end
% end
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
    % search radius = min(r, number of positive components)
    rF = min(r, npForward);
    rB = min(r, npBackward);
%     if npForward <=0 && npBackward <= 0
%         STOP = 3; % local maximum
%     elseif npForward <= 0 
%         STOP = 1; % all coordinates want to decrease
%     elseif npBackward <= 0
%         STOP = 2; % all coordinates want to increase
%     end
    fprintf('Number of positive forward and backward gradients: %d and %d \n', npForward, npBackward);

    % From Cezik 2008: adaptive step size to avoid non-concave regions
    if SL_beta_x < 0.5
        d = 3;
    elseif SL_beta_x < 0.65
        d = 2;
    else
        d = 1;
    end       
%         nChange = min(nChange, d*(npForward+npBackward)); % e.g. if there are only 2 positive gradient components, make nChange <= 2
    changeTolerancePlus = d * ones(nAgentGroups, nShifts);
    changeToleranceMinus = d * ones(nAgentGroups, nShifts);
    countChange = 0;
    while countChange < nChange
        if npForward <= 0 && npBackward <= 0
            break % max benefit <= 0
        end

        % Add agents
        if rF > 0
            randIndex = randi(rF);
            temp = IndForward(randIndex);
            [pRow, pCol] = ind2sub(size(forwardGradient), temp);

            x_trial(pRow, pCol) = x_trial(pRow, pCol) + 1;
            countChange = countChange + 1;
            changeTolerancePlus(pRow, pCol) = changeTolerancePlus(pRow, pCol) - 1;
            if changeTolerancePlus(pRow, pCol) == 0
                IndForward(randIndex) = []; % remove the chosen index
                npForward = npForward - 1;
            end
        end

        % Remove agents
        if rB > 0
            randIndex = randi(rB);
            temp = IndBackward(randIndex);
            [mRow, mCol] = ind2sub(size(backwardGradient), temp);

            x_trial(mRow, mCol) = x_trial(mRow, mCol) - 1;
            countChange = countChange + 1;
            changeToleranceMinus(mRow, mCol) = changeToleranceMinus(mRow, mCol) - 1;

            % if the component has been changed d times, or any further
            % reduction gives negative agents
            if changeToleranceMinus(mRow, mCol) == 0 || x_trial(mRow, mCol) == 0
                IndBackward(randIndex) = [];       
                npBackward = npBackward - 1;
            end
        end
    end
        
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
end % end local search for x

end

