function [ f, x_opt ] = OptimizeMultiSkill(  )
% Pot 2007
% Using  MultiSkill(x, runlength, seed, nCallTypes, nAgentGroups, varsigma, arrivalRates, R, serviceLevelMin)
clc;clear;
fileName = datestr(clock,'mm_dd_HH_MM_SS');
diary([pwd, '\outputs\', fileName, '.txt']);

% ------------------------ TAKEN FROM MultiSkill --------------------------
serviceLevelMin = 0.8; % service level requirement
meanST = 8; % mean service time
% Problem configuration (small call center)
nCallTypes =2;
nAgentGroups = 2; % types of agents who have different skills
arrivalRates = csvread('ArrivalRatesSmall.csv');
% Text file shifts contains all of the shifts being considered in the
% following format:
% col 1. Time shift starts, in minutes.
% col 2. Time of first break (15 min).
% col 3. Time of lunch break (30 min).
% col 4. Time of second break (15 min).
% col 5. Time shift ends
shifts = csvread('Shifts.csv');
% Defines the different agent groups and their preferences. In this case,
% there are two groups. the first one can only take calls of type 1 while
% the second one can take both types, yet it prefers type 2 calls. A value
% of 1000 means that calls cannot be taken. (Inf in problem statement, but
% 1000 is easier to handle in excel)
R = [1 1000; 2 1]; % nGroups x nCallTypes

% starting solution
nShifts = size(shifts, 1);
% x = zeros(nAgentGroups, nShifts);
% x(1,57) = 20;
% x(1,68) = 20;
% x(2,57) = 10;
% x(2,68) = 10;

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
M = 1e4; % arbitrary big number
epsilon = 0.1; % stopping precision for beta
max_iter = 1; % maximum number of consecutive fails in local search for x
nChangeMax = 32; %floor(nAgentGroups * nShifts); % maximum number of swaps in x when generating trial solution
golden_ratio = (1 + sqrt(5))/2;

% golden section search for n
fAll = [];
f = -M; % best solution so far
n_L = 0; % initial bounds for #agents
n_U = max(max(arrivalRates)) * meanST  / 15; % required by 15-min max workload required
count_n = 0;
n_sim = 0;
while n_U - n_L > 1
    n = ceil(n_U - (n_U - n_L) / golden_ratio); % new "center" point
    count_n = count_n + 1;
    fprintf('==================================================== \n');
    fprintf('Golden section iteration (%d) of n: n = %d. \n', count_n, n);
    fprintf('n_sim = %d. \n', n_sim);
    
    % starting by allocating agents uniformly among groups & shifts
    x = zeros(nAgentGroups, nShifts);
    for k = 1:n
        r = randi(nAgentGroups);
        c = randi(nShifts);
        x(r, c) = x(r, c) + 1;
    end
    x_ast = x;
    x_ast_ast = x;

    % bisection search for beta
    beta_L = 0; % initial bounds for beta
    beta_U = M;
    f_n = M; % best solution when number of agents = n
    count_beta = 0;
    while beta_U - beta_L > epsilon
        z = max(2, log(beta_U - beta_L));
        beta = beta_L + (beta_U - beta_L)/z; % new center point
        count_beta = count_beta + 1;
        fprintf('------------------------------------------------ \n');
        fprintf('Bisection iteration (%d) of beta: beta = %d. \n', count_beta, beta);
        
        % local search for all x such that sum(sum(x)) = n
%         [costAll, SLAll] = MultiSkill(x, beta, runlength, seed, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
        if count_beta == 1 % initialize
            [f_n_beta, SL, sd, forwardGradient, backwardGradient] = MultiSkillPickedCalls(x, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
            fprintf('Initialize: f_n_beta = %d +/- %d, SL = %d. \n', f_n_beta, sd, SL);
        end
            %     cost = mean(costAll);
    %     SL = mean(SLAll);
%         f_n_beta = beta * (SL - serviceLevelMin) - cost;
        n_sim = n_sim + runlength;
%         [forwardGradient, backwardGradient] = GradientTable(x, beta, nCallTypes, nCalls, meanST, lateCalls, lastCalls, R, shifts);
        backwardGradient(x==0) = NaN; % remove impossible reductions in x
        nChange = nChangeMax; 
        count_x = 0;
        count_failed_x = 0;
        stopped = 0;
        % stop when either consecutive fails exceed tolerance, or there's
        % no room for improvement
        while count_failed_x < max_iter
            x = x_ast; % restore last best solution
            count_x = count_x + 1;
            fprintf('******************************************** \n');
            fprintf('Local search iteration (%d) of x: nChange = %d. \n', count_x, nChange);
            
            countChange = 0;
            while countChange < nChange && ~stopped
                % Generate trial solution x: find maximum in forward and
                % backward gradients, and move one agent accordingly
                % Each group*shift combination will be changed by 1 at most
                [pMax, temp] = max(forwardGradient(:));
                [pRow, pCol] = ind2sub(size(forwardGradient), temp);
                forwardGradient(pRow, pCol) = -Inf;
                [mMax, temp] = max(backwardGradient(:));
                [mRow, mCol] = ind2sub(size(backwardGradient), temp);
                backwardGradient(mRow, mCol) = -Inf;
                if pMax + mMax < 0
                    stopped = 1; % net benefit < 0
                else
                    if x(mRow, mCol) > 0
                        x(mRow, mCol) = x(mRow, mCol) - 1;
                        x(pRow, pCol) = x(pRow, pCol) + 1;
                        countChange = countChange + 1;    
                    else
                        backwardGradient(mRow, mCol) = -Inf;
                    end
                end
            end
            
            [f_beta_s, SL_x, sd_x, forwardGradient, backwardGradient] = MultiSkillPickedCalls(x, beta, ceil(runlength/3), seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
%             [costAll, SLAll, forwardGradient, backwardGradient] = MultiSkill(x, runlength, seed, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
%             [forwardGradient, backwardGradient] = GradientTable(x, beta, nCallTypes, nCalls, meanST, lateCalls, lastCalls, R, shifts);
%             cost = mean(costAll);
%             SL_x = mean(SLAll);
            n_sim = n_sim + ceil(runlength/3);
%             f_beta_s = beta * (SL_x - serviceLevelMin) - cost;

            if f_beta_s > f_n_beta
                fprintf('x step: improve f_n_beta from %d to %d +/- %d. \n', f_n_beta, f_beta_s, sd_x);
                sd = sd_x;
                x_ast = x; % update last best solution
                f_n_beta = f_beta_s;
                SL = SL_x;
                count_failed_x = 0;
                nChange = nChangeMax;
            else
                count_failed_x = count_failed_x + 1;
                fprintf('x step: failed for %d times. \n', count_failed_x);
                if count_failed_x == max_iter && nChange > 1
                    nChange = floor(nChange/2);
                    count_failed_x = 0;
                end
            end
        end % end local search for x
        
        [f_n_beta, SL, sd, forwardGradient, backwardGradient] = MultiSkillPickedCalls(x_ast, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
        if f_n_beta < f_n
            fprintf('beta step: improve f_n from %d to %d +/- %d. ', f_n, f_n_beta, sd);
            x_ast_ast = x_ast;
            f_n = f_n_beta;
        else
            fprintf('beta step: failed. ');
        end
    
        % Bisection beta step: find beta s.t. SL = serviceLevelMin
        if SL < serviceLevelMin
            beta_L = beta;
        else
            beta_U = beta;
        end
        fprintf('beta_L = %d, beta_U = %d. \n', beta_L, beta_U);

    end % end bisection search for beta
    
    % Golden section n step:
    % always increase n if SL is not met
    if SL < serviceLevelMin
        n_L = n;
    end
    % update bounds and record best solution
    if f_n > f
        n_L = n;
        f = f_n;
        fAll = [fAll, [n_sim; n; beta; f]];
        x_opt = x_ast_ast;
        fprintf('n step: improve f from %d to %d +/- %d. n_L = %d, n_U = %d. \n', f, sd, f_n, n_L, n_U);
    else
        n_U = n;
        fprintf('n step: failed. n_L = %d, n_U = %d. \n', n_L, n_U);
    end
    
end

save([pwd, '\outputs\', fileName]);
[f_opt, SL_opt, sd_opt] = MultiSkillPickedCalls(x_opt, beta, 50, seed+1, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
n_sim = n_sim + 50;
% [cost, SL_x] = MultiSkill(x_opt, 50, seed+1, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
% f_opt = beta * (SL_x - serviceLevelMin) - cost;
fprintf('Number of simulations: %d (each with runlength = %d). \n', n_sim, runlength);
fprintf('Optimal cost: %d +/- %d. \n', f_opt, sd_opt);
fprintf('Optimal SL: %d. \n', SL_opt);
save([pwd, '\outputs\', fileName]);
diary; % diary off
        
end % end golden section search for n

