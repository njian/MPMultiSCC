clc;clear;
fileName = datestr(clock,'mm_dd_HH_MM_SS');
diary([pwd, '\outputs\', fileName, '.txt']);
fprintf('OptimizeMultiSkillNew \n');

%% Read Data
% Everything is in MINUTES.
% Problem configuration (small call center)
% nCallTypes = 2;
% nAgentGroups = 2; % types of agents who have different skills
% meanST = 1/8 * 60 * ones(nCallTypes, nAgentGroups); % service t in minutes
% arrivalRates = csvread('ArrivalRatesSmall.csv');
% arrivalRates = arrivalRates / 60; % convert into # per minute
% Defines the different agent groups and their preferences. In this case,
% there are two groups. the first one can only take calls of type 1 while
% the second one can take both types, yet it prefers type 2 calls. A value
% of 1000 means that calls cannot be taken. (Inf in problem statement, but
% 1000 is easier to handle in excel)
% R = [1 1000; 2 1]; % nGroups x nCallTypes

% Problem configuration (large call center)
% %%%%%%%%%%%%%%For larger call-center, comment out all previous %%%%%%%%%%%%%
% %%%%%%%%%%%%%%definitions except for reading shifts and uncom- %%%%%%%%%%%%%
% %%%%%%%%%%%%%%ment the following:                              %%%%%%%%%%%%%

% nCallTypes = 20;
% nAgentGroups = 35;
% R = csvread('Routing.csv');
% arrivalRates = csvread('ArrivalRatesLarge.csv');
% arrivalRates = arrivalRates / 60; % convert into # per minute
% meanST = 1/8 * 60 * ones(nCallTypes, nAgentGroups); % service t in minutes
% costByGroup = [];
% % 
% %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % -------------------------------------------------------------------------
% 
% % Text file shifts contains all of the shifts being considered in the
% % following format:
% % col 1. Time shift starts, in minutes.
% % col 2. Time of first break (15 min).
% % col 3. Time of lunch break (30 min).
% % col 4. Time of second break (15 min).
% % col 5. Time shift ends
% shifts = csvread('ShiftsAvramidis.csv');
% varsigma = 0.1;

% Problem configuration (Pot 2007, small call center)
nCallTypes = 3;
nAgentGroups = 6; % types of agents who have different skills
% Case ID 1
arrivalRates = [1.0, 1.5, 2.0]';
R = [1    1000 1000;
     1000 1    1000;
     2    2    1000;
     3    1000 1   ;
     1000 3    2   ;
     4    4    3   ]; % nGroups x nCallTypes
meanST = 1./[0.20 NaN NaN; NaN 0.18 NaN; 0.19 0.17 NaN; 0.19 NaN 0.16; NaN 0.17 0.16; 0.18 0.16 0.15]'; % minutes
costByGroup = [1.0, 1.0, 1.1, 1.1, 1.1, 1.2]';
shifts = [];

SLtime = 20/60; % 20 seconds
serviceLevelMin = 0.8;

%% Algorithm
% Parameters
seed = 1;
runlength = 1;
fprintf('seed = %d, runlength = %d. \n', seed, runlength);

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

% Generate two random streams, one for initial solution, one for local
% search.
[SearchStream] = RandStream.create('mrg32k3a', 'NumStreams', 1);
SearchStream.Substream = seed;
RandStream.setGlobalStream(SearchStream);
tic;

[ x_opt, beta, history_n, history_all, n_sim, n_beta ] = goldenSection_n( runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime, costByGroup);

% vary n, bisection beta
% [ ~, x_opt, ~, ~, beta, n_sim, history_all ] = bisection_varyN_beta(40, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, 0);
% fix n, bisection beta
% [~, x_opt, ~, ~, beta, n_sim, history_all] = bisection_beta(38, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, 1);
t = toc;

%% Final Outputs  
fprintf('************************************************************************ \n');
fprintf('************************************************************************ \n');
fprintf('************************************************************************ \n');
save([pwd, '\outputs\', fileName]);
if ~isnan(x_opt)
    [f_opt, SL_opt, sd_opt, ~, ~, avgAbandon, avgServerUtil] = MultiSkillPickedCallsPots(x_opt, beta, ceil(runlength/3*5), seed+1, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, SLtime, costByGroup);
%     [f_opt, SL_opt, sd_opt, ~, ~, avgAbandon, avgServerUtil] = MultiSkillPickedCalls(x_opt, beta, ceil(runlength/3*5), seed+1, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, SLtime);
    cost_opt = beta * (SL_opt - serviceLevelMin) - f_opt;
    fprintf('Number of simulated days: %d. Each solution evaluated with %d days. \n', n_sim, runlength);
    fprintf('Number of beta steps: %d. \n', n_beta);
    fprintf('Total CPU time: %d. \n', t);
    fprintf('beta = %.2f. \n', beta);
    fprintf('Optimal f: %.2f +/- %.2f. \n', f_opt, sd_opt);
    fprintf('Optimal cost: %.2f. \n', cost_opt);
    fprintf('Optimal SL: %.2f. \n', SL_opt);
    fprintf('Number of agents: %d. \n', sum(sum(x_opt)));
else
    fprintf('Error: all n steps failed.');
end
avgServerUtilByGroup = nanmean(avgServerUtil,2);
save([pwd, '\outputs\', fileName]);
diary; % diary off
beep;

%% Test Optimality
% fprintf('Testing local optimality using finite differences: \n');
% testOptimality
% save([pwd, '\outputs\', fileName, 'CheckOptimality']);