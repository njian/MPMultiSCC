function [ fAll ] = OptimizeMultiSkill2(  )
% Without golden search for n
clc;clear;
fileName = datestr(clock,'mm_dd_HH_MM_SS');
diary([pwd, '\outputs\', fileName, '.txt']);
fprintf('OptimizeMultiSkill2 \n');
% fprintf('Backward gradient = - Forward \n');


% ------------------------ TAKEN FROM MultiSkill --------------------------
% Text file shifts contains all of the shifts being considered in the
% following format:
% col 1. Time shift starts, in minutes.
% col 2. Time of first break (15 min).
% col 3. Time of lunch break (30 min).
% col 4. Time of second break (15 min).
% col 5. Time shift ends
shifts = csvread('Shifts.csv');
nShifts = size(shifts, 1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem configuration (small call center)
serviceLevelMin = 0.8; % service level requirement
meanST = 8; % mean service time
nCallTypes =2;
nAgentGroups = 2; % types of agents who have different skills
arrivalRates = csvread('ArrivalRatesSmall.csv');
% Defines the different agent groups and their preferences. In this case,
% there are two groups. the first one can only take calls of type 1 while
% the second one can take both types, yet it prefers type 2 calls. A value
% of 1000 means that calls cannot be taken. (Inf in problem statement, but
% 1000 is easier to handle in excel)
R = [1 1000; 2 1]; % nGroups x nCallTypes

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem configuration (large call center)
% %%%%%%%%%%%%%%For larger call-center, comment out all previous %%%%%%%%%%%%%
% %%%%%%%%%%%%%%definitions except for reading shifts and uncom- %%%%%%%%%%%%%
% %%%%%%%%%%%%%%ment the following:                              %%%%%%%%%%%%%
%
% nCallTypes = 20;
% nAgentGroups = 35;
% varsigma = 0.1;
% serviceLevelMin = 0.8;
% meanST = 8;
% patMean = 10;
% SLtime = 20;
% nDays = runlength;
%
% R = csvread('Routing.csv');
% arrivalRates = csvread('ArrivalRatesLarge.csv');
%
% schedule = x;
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Transform routing R into Route. Each row of Route lists the possible agent
%groups that can take the corresponding call type, in decreasing order of
%preference.
Rcopy = R; % nGroups x nCallTypes
Route = zeros(nCallTypes, nAgentGroups);
for j=1:nCallTypes
    [C,index] = min(Rcopy(:,j)); % most preferred group for each callType
    count = 1;
    while(C ~= 1000)
        Route(j,count) = index; % Route(callType, count-th) = count-th preferred corresponding group
        count = count+1;
        Rcopy(index,j)=1000; % mark Rcopy of that to be "visited"
        [C,index] = min(Rcopy(:,j)); % find the next one
    end
end
clear Rcopy; clear C;
% -------------------------------------------------------------------------

% search parameters
seed = 1;
runlength = 30;
M = 100; % arbitrary big number
epsilon = 0.1; % stopping precision for beta
nChangeMax = 1; %2^floor(log2(nAgentGroups*nShifts/4)); % maximum number of swaps in x when generating trial solution
r = 10; % local random search among this many largest gradient components
max_fail = 20; % maximum number of consecutive fails in local search for x

% generate two random streams, one for initial solution, one for local
% search
[InitializeStream, SearchStream] = RandStream.create('mrg32k3a', 'NumStreams', 2);
InitializeStream.Substream = seed;
SearchStream.Substream = seed;

n_sim = 0;
fAll = [];
f = 1e4; % best solution so far
cost = 1e4;
n = ceil(max(max(arrivalRates)) * meanST  / 15); % required by 15-min max workload required
% starting by allocating agents uniformly among groups & shifts
RandStream.setGlobalStream(InitializeStream);
x_trial = zeros(nAgentGroups, nShifts);
for k = 1:n
    row = randi(nAgentGroups);
    col = randi(nShifts);
    x_trial(row, col) = x_trial(row, col) + 1;
end
clear row; clear col;
x_ast = x_trial;
x_opt = x_trial;

%% Algorithm
% Outer layer: f = min_beta f_beta, with bisection
% Inner layer: f_beta = max_x f_beta_x, with gradient-like search
% f_beta_x = -cost + beta*(SL-MinServiceLevel) for any beta and x
% bisection search for beta
RandStream.setGlobalStream(SearchStream);
beta_L = 0; % initial bounds for beta
beta_U = M;
count_beta = 0;
count_x = 0;
while beta_U - beta_L > epsilon
    z = max(2, log2(beta_U - beta_L));
    beta = beta_L + (beta_U - beta_L)/z; % new center point
    count_beta = count_beta + 1;
    fprintf('Bisection iteration (%d) of beta: beta_L = %.2f, beta_U = %.2f. \n', count_beta, beta_L, beta_U);
    fprintf('beta = %.2f. \n', beta);
    
    if count_beta == 1 % initialize
        [f_beta_x, SL_beta_x, sd_beta_x, forwardGradient, backwardGradient] = MultiSkillPickedCalls(x_ast, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);       
%         [ forwardGradient, backwardGradient ] = GradientTableFiniteDiff( x_ast, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts );
        backwardGradient(x_ast==0) = -Inf;
        n_sim = n_sim + runlength;
        f_beta = f_beta_x;
        fprintf('Initialize: f_beta = %.2f +/- %.2f, SL = %.2f. \n', f_beta_x, sd_beta_x, SL_beta_x);
        fprintf('Starting n: n = %d. \n', n);
%     else
        % SL and cost from last best solution x_ast (with new beta)
%         f_beta_x = beta * (SL_beta - serviceLevelMin) - cost_beta;
%         f_beta = f_beta_x;
%         fprintf('Warm Start: f_beta = %.2f. \n', f_beta_x);
    end
    % "warm-start" local search from the optimal x from last beta
    x_trial = x_ast;
    
    nChange = nChangeMax; 
    count_failed_x = 0;
    % stop when either consecutive fails exceed tolerance, or there's
    % no room for improvement
%     sum(forwardGradient(:)>0) > 0 || sum(backwardGradient(:)>0) > 0
%     r = 5;sa
    STOP = 0;
    while count_failed_x < max_fail && STOP == 0
        count_x = count_x + 1;
        fprintf('------------------------------------------------ \n');
        fprintf('Local search iteration (%d) of x: nChange = %d. \n', count_x, nChange);
%         if nChange == 0
%            nChange 
%         end

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
        if npForward <= 0 && npBackward <= 0
            STOP = 1; % no room for improvement for this x_ast
        end
        fprintf('Number of positive forward and backward gradients: %d and %d \n', npForward, npBackward);
        % search radius = min(r, number of positive components)
        rF = min(r, npForward);
        rB = min(r, npBackward);

        % From Cezik 2008: adaptive step size to avoid non-concave regions
        if SL_beta_x < 0.5
            d = 3;
        elseif SL_beta_x < 0.65
            d = 2;
        else
            d = 1;
        end       
%         nChange = min(nChange, d*(npForward+npBackward)); % e.g. if there are only 2 positive gradient components, make nChange <= 2
        changeTolerance = d * ones(nAgentGroups, nShifts);
        
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
                changeTolerance(pRow, pCol) = changeTolerance(pRow, pCol) - 1;
                if changeTolerance(pRow, pCol) == 0
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
                changeTolerance(mRow, mCol) = changeTolerance(mRow, mCol) - 1;
                
                % if the component has been changed d times, or any further
                % reduction gives negative agents
                if changeTolerance(mRow, mCol) == 0 || x_trial(mRow, mCol) == 0
                    IndBackward(randIndex) = [];       
                    npBackward = npBackward - 1;
                end
            end
        end

        [f_beta_x, SL_beta_x, sd_x, forwardGradient_x, backwardGradient_x] = MultiSkillPickedCalls(x_trial, beta, ceil(runlength), seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
%         [ forwardGradient, backwardGradient ] = GradientTableFiniteDiff( x_trial, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts );
        fprintf('Trial solution obj = %.2f, SL = %.2f. \n', f_beta_x, SL_beta_x);
        
        n_sim = n_sim + runlength;

        if f_beta_x > f_beta
            fprintf('x step: improve f_beta from %.2f to %.2f +/- %.2f. \n', f_beta, f_beta_x, sd_x);
            fprintf('x step: n = %d. \n', sum(sum(x_trial)));
            count_failed_x = 0;
%             STOP = 0;
            % update last best solution for this beta
            x_ast = x_trial;
            f_beta = f_beta_x;
            sd_beta = sd_x;
            SL_beta = SL_beta_x;
            cost_beta = beta * (SL_beta - serviceLevelMin) - f_beta; 
            forwardGradient = forwardGradient_x;
            backwardGradient = backwardGradient_x;
            backwardGradient(x_ast==0) = -Inf; % remove impossible reductions in x
            fAll = [fAll, [n_sim; f_beta; sd_beta; cost_beta; SL_beta; beta; sum(sum(x_ast)); 1]];
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

    fprintf('************************************************************************ \n');
    if f_beta <= f %cost_beta <= cost && SL_beta >= serviceLevelMin
        fprintf('beta step: improve f from %.2f to %.2f +/- %.2f. \n', f, f_beta, sd_beta);
        % update the best solution for any beta
        x_opt = x_ast;
        fprintf('Number of agents: %d. \n', sum(sum(x_opt)));
        f = f_beta;
        SL = SL_beta;
        cost = beta * (SL - serviceLevelMin) - f; 
        fAll = [fAll, [n_sim; f; sd_beta; cost; SL; beta; sum(sum(x_opt)); 0]];
    else
        x_ast = x_opt;
        fprintf('beta step: failed: this solution %.2f > last best %.2f. \n', f_beta, f);
        fprintf('Number of agents remains at: %d. \n', sum(sum(x_opt)));
    end

    % Bisection beta step: find beta s.t. SL = serviceLevelMin
    if SL_beta < serviceLevelMin
        beta_L = beta;
    else
        beta_U = beta;
    end

    [f_beta, SL_beta, sd_beta, forwardGradient, backwardGradient] = MultiSkillPickedCalls(x_ast, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
    backwardGradient(x_ast==0) = -Inf;
    fprintf('End beta step: f_beta = %.2f, SL_beta = %.2f. \n', f_beta, SL_beta);
    n_sim = n_sim + runlength;
    
end % end bisection search for beta
    
fprintf('************************************************************************ \n');
fprintf('************************************************************************ \n');
fprintf('************************************************************************ \n');
save([pwd, '\outputs\', fileName]);
[f_opt, SL_opt, sd_opt] = MultiSkillPickedCalls(x_opt, beta_U, 50, seed+1, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
cost_opt = beta_L*(SL_opt-serviceLevelMin)-f_opt;
better = fAll(4,:)<=cost_opt;
meetSL = fAll(5,:)>=serviceLevelMin;
n_sim = n_sim + ceil(runlength/3*5);
fprintf('Number of simulations: %d (each with runlength = %d). \n', n_sim, runlength);
fprintf('Finishing beta_L = %.2f, beta_U = %.2f. \n', beta_L, beta_U);
fprintf('Optimal f: %.2f +/- %.2f. \n', f_opt, sd_opt);
fprintf('Optimal cost: %.2f. \n', beta_U * (SL_opt - serviceLevelMin) - f_opt);
fprintf('Optimal SL: %.2f. \n', SL_opt);
fprintf('Number of agents: %d. \n', sum(sum(x_opt)));
fprintf('Any x_ast on path better than x_opt? %d. \n', sum(better.*meetSL));
save([pwd, '\outputs\', fileName]);
diary; % diary off
        
end