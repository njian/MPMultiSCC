clc;clear;
fileName = datestr(clock,'mm_dd_HH_MM_SS');
diary([pwd, '\outputs\', fileName, '.txt']);
fprintf('OptimizeMultiSkillNew \n');

%% Read Data
% Problem configuration (small call center)
serviceLevelMin = 0.8; % service level requirement
meanST = 8; % mean service time
% Problem configuration (small call center)
nCallTypes = 2;
nAgentGroups = 2; % types of agents who have different skills
arrivalRates = csvread('ArrivalRatesSmall.csv');
% Defines the different agent groups and their preferences. In this case,
% there are two groups. the first one can only take calls of type 1 while
% the second one can take both types, yet it prefers type 2 calls. A value
% of 1000 means that calls cannot be taken. (Inf in problem statement, but
% 1000 is easier to handle in excel)
R = [1 1000; 2 1]; % nGroups x nCallTypes

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
% 
% R = csvread('Routing.csv');
% arrivalRates = csvread('ArrivalRatesLarge.csv');

%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------

% Text file shifts contains all of the shifts being considered in the
% following format:
% col 1. Time shift starts, in minutes.
% col 2. Time of first break (15 min).
% col 3. Time of lunch break (30 min).
% col 4. Time of second break (15 min).
% col 5. Time shift ends
shifts = csvread('Shifts.csv');

%% Algorithm
% Parameters
seed = 1;
runlength = 1;

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
[SearchStream] = RandStream.create('mrg32k3a', 'NumStreams', 2);
SearchStream.Substream = seed;
RandStream.setGlobalStream(SearchStream);

[ x_opt, beta, fAll, n_sim ] = goldenSection_n( runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts);
% [ ~, x_opt, ~, ~, beta, n_sim ] = bisection_beta(86, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts, 0);

%% Final Outputs  
fprintf('************************************************************************ \n');
fprintf('************************************************************************ \n');
fprintf('************************************************************************ \n');
save([pwd, '\outputs\', fileName]);
if ~isnan(x_opt)
    [f_opt, SL_opt, sd_opt] = MultiSkillPickedCalls(x_opt, beta, ceil(runlength/3*5), seed+1, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
    n_sim = n_sim + ceil(runlength/3*5);
    fprintf('Number of simulations: %d (each with runlength = %d). \n', n_sim, runlength);
    fprintf('Finishing beta_U = %.2f. \n', beta);
    fprintf('Optimal f: %.2f +/- %.2f. \n', f_opt, sd_opt);
    fprintf('Optimal cost: %.2f. \n', beta * (SL_opt - serviceLevelMin) - f_opt);
    fprintf('Optimal SL: %.2f. \n', SL_opt);
    fprintf('Number of agents: %d. \n', sum(sum(x_opt)));
else
    fprintf('Error: all n steps failed.');
end
save([pwd, '\outputs\', fileName]);
diary; % diary off