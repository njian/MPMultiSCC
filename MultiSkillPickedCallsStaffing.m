function [f, SL, stdev, forwardGradient, backwardGradient, avgAbandon, avgServerUtil, lateCalls, lastCalls1, lastCalls2, lastCalls3, totalCalls] = MultiSkillPickedCallsStaffing(x, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, SLtime, shifts, costByGroup)
% x is a matrix (nGroups x 1) containing the number of agents in group i 
% beta is the Lagrangian multiplier for the SL constraint
% runlength is the number of days of simulated time to simulate
% seed is the index of the substreams to use (integer >= 1)
% serviceLevelMin is the minimal overall service level costraint
% nCallTypes is the number of types of calls or skills
% nAgentGroups is the number of agent groups or skill-sets
% arrivalRates is by nCallTypes times (15 minute) periods, in # per minute
% each period.
% meanST is nCallTypes by nAgentGroups mean service time, in MINUTES.
% R is nAgentGroups by nCallTypes, representing types of calls feasible for
% each agent group. Smaller means more preferred and 1000 means infeasible.
% Route is nCallTypes by nAgentGroups, representing the most to least
% preferred agent groups for each call type.

% RETURNS: Mean and confidence interval of the lowest (per period) % of calls
% answered in less than 20 seconds (should be >0.8) for each call type
% and total daily cost. To modify this output look at nCall and
% nLessThan20.

%   ***************************************
%   *** Code written by German Gutierrez***
%   ***         gg92@cornell.edu        ***
%   *** Updated by Danielle Lertola to  ***
%   *** use standard calling and random ***
%   *** number streams - 6/8/2012       ***
%   ***      Edited by Bryan Chong      ***
%   ***  bhc34@cornell.edu Nov 3rd 2014 ***
%   *** Modified by Nanjing Jian for    ***
%   *** optimization purpose            ***
%   ***  nj227@cornell.edu Apr 2017     ***
%   ***************************************

shifts = [];
SLAll = zeros(runlength, 1);
forwardGradient = zeros(nAgentGroups, 1);
backwardGradient = zeros(nAgentGroups, 1);
avgAbandon = zeros(nCallTypes, 1);
agentBusyTime = zeros(nAgentGroups, 1); 
t = 0;
%MULTI-SKILL CALL CENTER.
% Uses discrete event simulation with 4 events:
% 1. Call Arrives
% 2. Call is completed
% 3. Agent Starts/returns from break (sub code: 1 - starts, 2-returns from
%    first 15-min break, 3-returns from 30-min, 4-returns from 2nd 15min.)
% 4. Agent goes on break/leaves (sub code: 1- first 15min break,
%    2-30min break, 3- second 15-min and 4-completes shift.)
%
% Throughout the entire simulation, each agent and call is tracked
% independently to ensure that waiting times and breaks are treated
% correctly.
nDays = runlength;

if (sum(sum(x < 0)) > 0|| (sum(size(x) ~= [nAgentGroups, 1])>0) || (runlength <= 0) || (seed <= 0) || (round(seed) ~= seed))
    fprintf('x should be %d * 74 matrix with nonnegative components\nrunlength should be positive and real,\nseed should be a positive integer\n', nCallTypes);
    fn = NaN;
    FnVar = NaN;
else % main simulation
        
    % the number of agents in each group being hired.
    schedule = x; % nGroups * 1, number of agents to hire
   
    maxTime = 9600; % 160 hrs
    
    % GENERATE RANDOM NUMBER STREAMS
    % Generate new streams for call arrivals, call
    [ArrivalStream, PatientStream, ServiceStream, GradientStream] = RandStream.create('mrg32k3a', 'NumStreams', 4);
    % Set the substream to the "seed"
    ArrivalStream.Substream = seed;
    PatientStream.Substream = seed;
    ServiceStream.Substream = seed;
    GradientStream.Substream = seed;
    rng(seed);
    
    OldStream = RandStream.getGlobalStream(); % Temporarily store old stream
    
    for i=1:nDays
        % Calls that are picked up late (>20s) or expired without being picked up
        % each column: [call arrival time, call type, call expiry time, pick up time or Inf, service time or NaN]'
        lateCalls = [];
        % Calls that are picked up by an agent who is the last available one in the group since then
        % each column: [call arrival time, call type, agent group, agent shift, pick up time]'
        lastCalls1 = zeros(nAgentGroups, 1); lastCalls2 = zeros(nAgentGroups, 1); lastCalls3 = zeros(nAgentGroups,1);
        lastAgent1 = -1 * ones(nAgentGroups, 1); % the index of the agent in each group * shift
        lastAgent2 = -1 * ones(nAgentGroups, 1); % the index of the agent in each group * shift
        lastAgent3 = -1 * ones(nAgentGroups, 1); % the index of the agent in each group * shift

        %Events will track all future events:
        % call arrivals have (type 1, time,call type and Time+patience time)
        % call exits have (type 2, time, agent#, 0)
        % agent starts (type 3, time, agent #, starting subcode)
        % agent leaves (type 4,time, agent #, leave subcode).
        
        Events = zeros(4,nCallTypes); % nCallTypes because at least two arrivals would be scheduled. This will change size constantly after that.
        callQueue = zeros(3,1); % (time call entered queue; call type; call expiry time), equivalent to Events(2:4,i) for an arrival event (type 1).
        nPerGroup = sum(schedule,2); % nAgentGroups x 1 col vector with total number of agents in each group
        
        %Contain the number of calls that arrive/are answered in less than 20 seconds
        % for each period and type.
        nCalls= zeros(nCallTypes,1);
        nLessThan20 = zeros(nCallTypes,1);
        
        %Generate first calls via thinning method for non-homog poisson process
        for j=1:nCallTypes
            lambda = arrivalRates(j);
            nextTime = exprnd(1/lambda);
            Events(:,j) = [1; nextTime; j; 0]; % Generate patience times later
        end
        
        %Generate patience times
        for j=1:nCallTypes
            Events(4,j) = Inf;
        end
        
        %Transform schedule into Agent matrix and generate agent start events.
        % (Agent Group, Agent Shift, (1 for avail, 0 for not avail), time of
        % next availability)
        Agents = zeros(sum(sum(schedule)), 4); % total agents * 4
        count = 1;
        nAvailable = zeros(nAgentGroups,1);
        s = size(schedule);
        for l=1:s(1)         % Agent Group
            for j=1:s(2)       % Agent Shift
                for k=1:schedule(l,j)
                    Agents(count,:)= [l, j, 0, inf];
                    %Type 3: Agent starts: type 1
                    Events = [Events [ 3; 0; count; 1]]; % record all agent starts into events
                    count = count + 1;
                end
            end
        end
        
        [~, nextEvent] = min(Events(2,:)); % index of the soonest event
        time = Events(2,nextEvent);
        done = 0;
        prevTime = 0;
        while( done == 0 )
            %Check if all events have been dealt with, time >= 540
            if (prevTime >= maxTime)
                done = 1;
            end
            prevTime = time; % last executed event's time
            
            %Only accept calls until 540, all other events continue to occur
            %after 540 since calls already in the system must be serviced.
            if( Events(1,nextEvent) == 1 )
                if ( time <= maxTime )
                    % CALL ARRIVES
                    callType = Events(3,nextEvent);
                    nCalls(callType) = nCalls(callType) + 1; % number of calls by type and shift
                    count = 1;
                    nextGroup = Route(callType,count);
                    callTaken = 0;
                    while nextGroup ~= 0 && count <= nAgentGroups % go through each feasible agent group
                        if nAvailable(nextGroup) > 0 % if there's available agent in that group
                            %Set stream for generating service time of agent
                            RandStream.setGlobalStream(ServiceStream);
                            %find agent, takes the call
                            for j = 1+sum(nPerGroup(1:nextGroup-1)):sum(nPerGroup(1:nextGroup)) % all agents in nextGroup
                                if( Agents(j,4) <= time && Agents(j,3) == 1 ) % next available time is now and agent is available
                                    sTime = exprnd(meanST(callType, nextGroup)); % generate service time
                                    agentBusyTime(Agents(j,1), Agents(j,2)) = agentBusyTime(Agents(j,1), Agents(j,2)) + sTime;
                                    Agents(j,4) = time + sTime; % update availability
                                    Agents(j,3) = 0;
                                    nAvailable(nextGroup) = nAvailable(nextGroup) - 1;
                                    if nAvailable(nextGroup) == 2 && lastAgent3(nextGroup, Agents(j,2)) == -1 % when no agent is labeled and this is the third to last
                                        lastAgent3(nextGroup, Agents(j,2)) = j; % label this agent
                                        lastCalls3(nextGroup, Agents(j,2)) = lastCalls3(nextGroup, Agents(j,2)) + 1; % record this call
                                    elseif j == lastAgent3(Agents(j,1),Agents(j,2)) % when this is labeled as the third to last agent
                                        lastCalls3(nextGroup, Agents(j,2)) = lastCalls3(nextGroup, Agents(j,2)) + 1; % record this call
                                    end
                                    if nAvailable(nextGroup) == 1 && lastAgent2(nextGroup, Agents(j,2)) == -1 % when no agent is labeled and this is the second to last
                                        lastAgent2(nextGroup, Agents(j,2)) = j; % label this agent
                                        lastCalls2(nextGroup, Agents(j,2)) = lastCalls2(nextGroup, Agents(j,2)) + 1; % record this call
                                    elseif j == lastAgent2(Agents(j,1),Agents(j,2)) % when this is labeled as the second to last agent
                                        lastCalls2(nextGroup, Agents(j,2)) = lastCalls2(nextGroup, Agents(j,2)) + 1; % record this call
                                    end
                                    if nAvailable(nextGroup) == 0 && lastAgent1(nextGroup, Agents(j,2)) == -1 % when no agent is labeled and this is the last
                                        lastAgent1(nextGroup, Agents(j,2)) = j; % label this agent
                                        lastCalls1(nextGroup, Agents(j,2)) = lastCalls1(nextGroup, Agents(j,2)) + 1; % record this call
                                    elseif j == lastAgent1(Agents(j,1),Agents(j,2)) % when this is labeled as the last agent
                                        lastCalls1(nextGroup, Agents(j,2)) = lastCalls1(nextGroup, Agents(j,2)) + 1; % record this call
                                    end
                                    Events = [Events [2; time+sTime; j; 0]]; % schedule call exit
                                    nextGroup = 0;
                                    callTaken = 1;
                                    %Call answered immediately: success
                                    nLessThan20(callType) = nLessThan20(callType) + 1;
                                    break
                                end
                            end
                        else
                            count= count+1;
                            if (count <= nAgentGroups)
                                nextGroup = Route(callType,count);
                            end
                        end
                    end
                    %No agent immediately available. Enters call queue
                    if callTaken == 0
                        callQueue =[callQueue Events(2:4,nextEvent)];
                    end
                    
                    RandStream.setGlobalStream(ArrivalStream);
                    %Generate next call arrival of this type.
                    success = 0;
                    lambda = arrivalRates(callType);
                    arrTime = time + exprnd(1/lambda);
                    
                    %Generate Patience time
                    RandStream.setGlobalStream(PatientStream);
                    Events = [Events [1; arrTime; callType; Inf]];
                end
                Events(:,nextEvent) = [];
                [~, nextEvent] = min(Events(2,:));
                time = Events(2,nextEvent);
            elseif(Events(1,nextEvent) == 2)
                % CALL IS COMPLETED:
                % If agent went on break, time of next availability, i.e.
                % Agents(i,4)  is not equal to current time.
                a = Events(3,nextEvent); % agent number who handled the call
                if Agents(a,4) == time % next availability is now
                    %Agent becomes available: Search for a call in queue and,
                    %if none, set Agent as idle.
                    j=2; callTaken = 0;
                    nCallQueue = size(callQueue,2);
                    while j <= nCallQueue && callTaken == 0
                            if (R(Agents(a,1),callQueue(2,j)) ~= 1000) % take first feasible call in queue
                                %Agent can take this call and abandonment time has
                                %not expired.
                                %Generate Service Time
                                RandStream.setGlobalStream(ServiceStream);
                                sTime = exprnd(meanST(callQueue(2,j),Agents(a,1)));
                                agentBusyTime(Agents(a,1), Agents(a,2)) = agentBusyTime(Agents(a,1), Agents(a,2)) + sTime;
                                Agents(a,4) = time + sTime;
                                Agents(a,3) = 0;
                                Events = [Events [2; Agents(a,4); a; 0]]; % schedule call exit
                                
                                %Count call as successful if answered withing
                                %SL time constraint.
                                if( time - callQueue(1,j) <= SLtime )
                                    nLessThan20(callQueue(2,j)) = nLessThan20(callQueue(2,j)) + 1;
                                else
                                    % call cannot be answered in time
                                    lateCalls = [lateCalls, [callQueue(:,j); time; sTime]];
                                end
                                
                                callQueue(:,j) = [];
                                nCallQueue = nCallQueue - 1;
                                callTaken = 1;
                            else
                                j = j+1;
                            end
%                         end
                    end
                    
                    if callTaken == 0
                        %Agent did not take any calls
                        Agents(a,3) = 1; % agent available
                        nAvailable(Agents(a,1)) = nAvailable(Agents(a,1)) + 1;
                    end
                end
                
                Events(:,nextEvent) = [];
                [~, nextEvent] = min(Events(2,:));
                time = Events(2,nextEvent);
                             
            elseif(Events(1,nextEvent) == 3)
                % AGENT STARTS WORKING
                %To do routing correctly, all agents starting at the same time
                %must be added to the available agents before assigning calls
                %to them.
                % First add all agents starting at the same time.     
                while Events(2,nextEvent) == time && Events(1,nextEvent) == 3 % agent arrival
                    a = Events(3,nextEvent); % agent group number
                    Agents(a,3) = 1; % set availability = 1
                    Agents(a,4) = Events(2,nextEvent); % next avail = call exit time
                    nAvailable(Agents(a,1)) = 1 + nAvailable(Agents(a,1)); % record number of available agent in that group
                    Events(:,nextEvent) = []; % clear this executed event
                    [~, nextEvent] = min(Events(2,:)); % find next soonest event in FEL
                end  
                time = Events(2,nextEvent);
            end
        end
        
        totalCalls = sum(sum(nCalls));
        SLAll(i,1) = sum(nLessThan20) / totalCalls;
        % Calculate gradient table
        RandStream.setGlobalStream(GradientStream)
        tic;
        [fg, bg] = GradientTableStaffing( beta, nAgentGroups, totalCalls, meanST, lateCalls, lastCalls1, R, costByGroup );
        t = t + toc;
        forwardGradient = forwardGradient + fg;
        backwardGradient = backwardGradient + bg;
        avgAbandon = avgAbandon ./ totalCalls;
        
        %Set arrival stream for next day
        RandStream.setGlobalStream(ArrivalStream);
    end
    %Return old Random Number Stream
    RandStream.setGlobalStream(OldStream);  
    
    cost = costByGroup' * x;
    obj = beta * (SLAll - serviceLevelMin) - cost;
    f = mean(obj);
    stdev = std(obj);
    SL = mean(SLAll);
    forwardGradient = forwardGradient / runlength;
    backwardGradient = backwardGradient / runlength;
    
    avgAbandon = avgAbandon / runlength;
    agentBusyTime = agentBusyTime ./ x / runlength;
    avgServerUtil = agentBusyTime ./ maxTime;
    fprintf('time to evaluate pseudo-gradient: %d.\n', t);
end
