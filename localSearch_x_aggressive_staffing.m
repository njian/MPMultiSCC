function [f_beta, x_ast, SL_beta, sd_beta, n_sim] = localSearch_x_aggressive_staffing(n, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime, costByGroup, x_ast)
% Given n and beta, local search for x.
if ~isempty(shifts)
    nShifts = size(shifts, 1);
else
    nShifts = 1;
end
 
% Search Parameters
d = 1; % max number of change to each coordinate of x
nChangeMax = 1; %min(32, ceil(nAgentGroups*nShifts*0.5*d)); % maximum number of swaps in x when generating trial solution 2^floor(log2(nAgentGroups*nShifts))
r = 1; %5; % local random search among this many largest gradient components
max_fail = 5; % maximum number of consecutive fails in local search for x
nRestarts = 3; % maximum number of restarts of the local search

count_x = 0;
nChange0 = nChangeMax;
count_failed_x = 0;
n_sim = 0;
triedPairs = [];
 
% Initialize
if isnan(x_ast)
    x_trial = evenly_spread( n, nAgentGroups, nShifts );
else
    x_trial = x_ast;
end
[f_beta_x, SL_beta_x, sd_x, forwardGradient_x, backwardGradient_x] = MultiSkillPickedCallsPots(x_trial, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, SLtime, costByGroup);
% [f_beta_x, SL_beta_x, sd_x, forwardGradient_x, backwardGradient_x] = MultiSkillPickedCalls(x_trial, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime);
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
while count_failed_x < max_fail && nChange0 > 0 && STOP == 0
    count_x = count_x + 1;
    fprintf('--------------------------------------------------------------------- \n');
    fprintf('Local search iteration (%d) of x. \n', count_x);
    fprintf('n = %d, beta = %.2f. \n', sum(sum(x_trial)), beta);
 
    if count_failed_x == 0
        % Sort the gradients if new gradients are obtained
        [SortedForward, IndForward] = sort(forwardGradient(:),'descend');
        [SortedBackward, IndBackward] = sort(backwardGradient(:),'descend');
    end
%     SortedForward = SortedForwardSave;
%     SortedBackward = SortedBackwardSave;
%     IndForward = IndForwardSave;
%     IndBackward = IndBackwardSave;
    % the positive components in forwardGradient and non-negative
    % components in backwardGradient indicate room for improvement
    npForward = sum(SortedForward>=0);
    npBackward = sum(SortedBackward>=0);
%     npForward0 = npForward;
%     npBackward0 = npBackward;
    maxRepetitiveChange = npForward0*npBackward0 + npForward0*(npBackward0==0) + npBackward0*(npForward0==0);
    nChange = min(nChange0, maxRepetitiveChange*d);
    if isempty(triedPairs)
        triedPairs = zeros(nChange,2);
    end
    if npForward <=0 && npBackward <= 0 
        STOP = 1; % local maximum
    end
    fprintf('Number of positive forward and backward gradients: %d and %d \n', npForward, npBackward);
    fprintf('nChange = %d. \n', nChange);
    
    % search radius <= number of positive backward/forward gradient
    % components so that there exists a forward/backward gradient component
    % to make the sum > 0.
    rB = min(r, sum(SortedBackward > -SortedForward(1)));
    rF = min(r, sum(SortedForward > -SortedBackward(1)));
    
    changeTolerancePlus = d * ones(nAgentGroups, nShifts);
    changeToleranceMinus = d * ones(nAgentGroups, nShifts);
    countChange = 0;
    countRepetitiveChange = 0;
    changedIndex = zeros(nChange,2);
    while countChange < nChange && rF > 0 && rB > 0 && countRepetitiveChange <= maxRepetitiveChange && STOP == 0 %&& (npForward > 0 || npBackward > 0)
        % Generate random group and shift to add an agent
        randIndexF = randi(rF);
        indF = IndForward(randIndexF);
        pCol = ceil(indF/nAgentGroups);
        pRow = indF - (pCol-1)*nAgentGroups;
        % Generate random group and shift to remove an agent, so that the
        % total change in obj is > 0
        max_rB = sum(SortedBackward > - forwardGradient(pRow, pCol));
        if max_rB > 0
            randIndexB = randi(max_rB);
            indB = IndBackward(randIndexB);
            mCol = ceil(indB/nAgentGroups);
            mRow = indB - (mCol-1)*nAgentGroups;            
            
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
            changedIndex(countChange,:) = [randIndexF, randIndexB];

            % update search radius
            if numel(SortedForward) > 0
                rB = min(r, sum(SortedBackward > -SortedForward(1)));
            else
                rB = 0;
            end
            if numel(SortedBackward) > 0
                rF = min(r, sum(SortedForward > -SortedBackward(1)));
            else
                rF = 0;
            end
        end
        
        % Check if the generated solution has been tried before
        if ~(countChange < nChange && rF > 0 && rB > 0 && countRepetitiveChange <= maxRepetitiveChange && STOP == 0)
            c = 1; match = 0; tried = 0;
            while tried == 0 && c <= size(triedPairs,3)
                for i = 1:size(changedIndex,1)
                    if changedIndex(i,1) > 0
                        match =  match + ismember(changedIndex(i,:), triedPairs(:,:,c), 'rows');
                    end
                end
                if match > 0 && match == sum(sum(triedPairs(:,:,c)>0))/2
                    tried = 1;
                end
                c = c + 1;
            end

            if tried == 0 
                % record pair
                try
                    triedPairs(:,:,count_failed_x+1) = changedIndex;
                catch
                    triedPairs
                end
                break; % stop generating trial solution
            else
                countRepetitiveChange = countRepetitiveChange + 1;
                x_trial = x_ast;
                countChange = 0; % go back and generate more trial solutions
            end    
        end
    end

    if countChange == 0  % local optimum by gradient, or by exhausting all tries
       fprintf('x step: no trial solution can be generated: rB = %d, rF = %d. \n', rB, rF); 
       STOP = 4;
    end
     
    if STOP == 1
        if nRestarts > 0
            fprintf('x step: Seems like a local minimum. nRestarts = %d. Restarting... \n', nRestarts);
            fprintf('--------------------------------------------------------------------- \n');
            n = sum(sum(x_trial));
            x_trial = ones(nAgentGroups, nShifts) * floor(n/(nAgentGroups*nShifts));   
            for k = 1:(n-sum(sum(x_trial)))
                row = randi(nAgentGroups);
                col = randi(nShifts);
                x_trial(row, col) = x_trial(row, col) + 1;
            end
            STOP = 0;
            count_failed_x = 0;
            nRestarts = nRestarts - 1;
        end
    end
    
    if STOP == 0
        fprintf('Evaluating trial solution.\n');
        [f_beta_x, SL_beta_x, sd_x, forwardGradient_x, backwardGradient_x] = MultiSkillPickedCallsPots(x_trial, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, SLtime, costByGroup);
%         [f_beta_x, SL_beta_x, sd_x, forwardGradient_x, backwardGradient_x] = MultiSkillPickedCalls(x_trial, beta, ceil(runlength), seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime);
        fprintf('Trial solution obj = %.2f, SL = %.2f. \n', f_beta_x, SL_beta_x);  
        n_sim = n_sim + runlength;
 
        if f_beta_x > f_beta
            fprintf('x step: IMPROVE f_beta from %.2f to %.2f +/- %.2f. \n', f_beta, f_beta_x, sd_x);
            count_failed_x = 0;
 
            % update last best solution for this beta
            x_ast = x_trial;
            f_beta = f_beta_x;
            sd_beta = sd_x;
            SL_beta = SL_beta_x;
            forwardGradient = forwardGradient_x;
            backwardGradient = backwardGradient_x;
            backwardGradient(x_ast==0) = -Inf; % remove impossible reductions in x
            triedPairs = [];
        else
            x_trial = x_ast;
            count_failed_x = count_failed_x + 1;
            fprintf('x step: failed for %d times: this solution %.2f < last best %.2f. \n', count_failed_x, f_beta_x, f_beta);
            % update stepsize
            if count_failed_x == max_fail && nChange0 > 1
                nChange0 = floor(nChange0/2);
                count_failed_x = 0;
                triedPairs = [];
            end
 
        end
        
    elseif STOP == 1
        % First parse found no room for improvement. Need to
        % increase/decrease n.
        fprintf('x step: no room for improvement. Seems like a local optimal. \n');
    end
end % end local search for x
 
end