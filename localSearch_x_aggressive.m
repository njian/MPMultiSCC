function [f_beta, x_ast, SL_beta, sd_beta, n_sim] = localSearch_x_aggressive(n, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts)
% Given n and beta, local search for x.
nShifts = size(shifts, 1);
 
% Search Parameters
nChangeMax = 32; %min(32, 2^floor(log2(nAgentGroups*nShifts/10))); % maximum number of swaps in x when generating trial solution
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
 
    if count_failed_x == 0
        % Sort the gradients if new gradients are obtained
        [SortedForwardSave, IndForwardSave] = sort(forwardGradient(:),'descend');
        [SortedBackwardSave, IndBackwardSave] = sort(backwardGradient(:),'descend');
    end
    SortedForward = SortedForwardSave;
    IndForward = IndForwardSave;
    SortedBackward = SortedBackwardSave;
    IndBackward = IndBackwardSave;
    % the positive components in forwardGradient and non-negative
    % components in backwardGradient indicate room for improvement
    npForward = sum(SortedForward>=0);
    npBackward = sum(SortedBackward>=0);
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
    
    % search radius <= number of positive backward/forward gradient
    % components so that there exists a forward/backward gradient component
    % to make the sum > 0.
    rB = min(r, sum(SortedBackward > -SortedForward(1)));
    rF = min(r, sum(SortedForward > -SortedBackward(1)));
    
    changeTolerancePlus = d * ones(nAgentGroups, nShifts);
    changeToleranceMinus = d * ones(nAgentGroups, nShifts);
    countChange = 0;
    q = 0;
    while countChange < nChange && rF > 0 && rB > 0 && STOP == 0 %&& (npForward > 0 || npBackward > 0)
        q = q + 1;
        % Generate random group and shift to add an agent
        randIndexF = randi(rF);
        [pRow, pCol] = ind2sub([nAgentGroups, nShifts], IndForward(randIndexF));
        % Generate random group and shift to remove an agent, so that the
        % total change in obj is > 0
        max_rB = sum(SortedBackward > - forwardGradient(pRow, pCol));
        if max_rB > 0
            randIndexB = randi(max_rB);
            [mRow, mCol] = ind2sub([nAgentGroups, nShifts], IndBackward(randIndexB));
            
            % Execute swap
            x_trial(pRow, pCol) = x_trial(pRow, pCol) + 1;
            changeTolerancePlus(pRow, pCol) = changeTolerancePlus(pRow, pCol) - 1;
            if changeTolerancePlus(pRow, pCol) <= 0
                IndForward(randIndexF) = []; % remove the chosen index
                SortedForward(randIndexF) = [];
            end

            x_trial(mRow, mCol) = x_trial(mRow, mCol) - 1; 
            changeToleranceMinus(mRow, mCol) = changeToleranceMinus(mRow, mCol) - 1;
            if (changeToleranceMinus(mRow, mCol) <= 0 || x_trial(mRow, mCol) <= 0)
                IndBackward(randIndexB) = [];
                SortedBackward(randIndexB) = [];
            end

            countChange = countChange + 1;
        end
        
        if q > 10
           q 
        end
        
        % search radius
        rB = min(r, sum(SortedBackward > -SortedForward(1)));
        rF = min(r, sum(SortedForward > -SortedBackward(1)));       
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
        fprintf('x step: no room for improvement. Seems like a local optimal. \n');
    end
end % end local search for x
 
end