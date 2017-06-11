function [f, SL, stdev, forwardGradient, backwardGradient] = MultiSkillPickedCalls(x, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, meanST, R, Route, shifts)
% function [fn, FnVar, FnGrad, FnGradCov, constraint, ConstraintCov, ConstraintGrad, ConstraintGradCov] = MultiSkill(x, runlength, seed, other);
% x is a matrix (nGroups x nShifts) containing the number of agents in
% group i working shift j
% beta is the Lagrangian multiplier for the SL constraint
% runlength is the number of days of simulated time to simulate
% seed is the index of the substreams to use (integer >= 1)
% serviceLevelMin is the minimal overall service level costraint
% nCallTypes is the number of types of calls or skills
% nAgentGroups is the number of agent groups or skill-sets
% arrivalRates is by nCallTypes times (15 minute) periods, in # per HOUR in
% each period.
% meanST is nCallTypes by nAgentGroups mean service time, in MINUTES.
% R is nAgentGroups by nCallTypes, representing types of calls feasible for
% each agent group. Smaller means more preferred and 1000 means infeasible.
% Route is nCallTypes by nAgentGroups, representing the most to least
% preferred agent groups for each call type.
% Shifts is the schedule of shifts by starting time, time of first break,
% time of lunch, time of second break, and ending time.

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
%   ***  nj227@cornell.edu Apr ??? 2017 ***
%   ***************************************

% Calls that are picked up late (>20s) or expired without being picked up
% each column: [call arrival time, call type, call expiry time, pick up time or Inf, service time or NaN]'
lateCalls = [];
% Calls that are picked up by an agent who is the last available one in the group since then
% each column: [call arrival time, call type, agent group, agent shift, pick up time]'
lastCalls = [];
lastAgent = -1 * ones(nAgentGroups, size(shifts, 1)); % the index of the agent in each group * shift
% pickedCalls(i,j) = callType if in minute i agent j picked up a call of
% this callType
nShifts = size(shifts, 1);
SLAll = zeros(runlength, 1);
forwardGradient = zeros(nAgentGroups, nShifts);
backwardGradient = zeros(nAgentGroups, nShifts);

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

SLtime = 20; % service level required time
patMean = 10; % patience exp with this mean
nDays = runlength;

if (sum(sum(x < 0)) > 0|| (sum(size(x) ~= [nAgentGroups, nShifts])>0) || (runlength <= 0) || (seed <= 0) || (round(seed) ~= seed))
    fprintf('x should be %d * 74 matrix with nonnegative components\nrunlength should be positive and real,\nseed should be a positive integer\n', nCallTypes);
    fn = NaN;
    FnVar = NaN;
else % main simulation
        
    % I x (number of shifts) matrix with the number of agents
    % in each group and schedule being hired.
    schedule = x; % nGroups * nShift, number of agents to hire
   
    maxTime = 540; % 9 hours
%     pickedCalls = sparse(maxTime, sum(sum(x)));
    
    % GENERATE RANDOM NUMBER STREAMS
    % Generate new streams for call arrivals, call
    [ArrivalStream, PatientStream, ServiceStream] = RandStream.create('mrg32k3a', 'NumStreams', 3);
    % Set the substream to the "seed"
    ArrivalStream.Substream = seed;
    PatientStream.Substream = seed;
    ServiceStream.Substream = seed;
    
    OldStream = RandStream.setGlobalStream(ArrivalStream); % Temporarily store old stream
    
    for i=1:nDays
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
        nCalls= zeros(nCallTypes, 36);
        nLessThan20 = zeros(nCallTypes,36);
        
        %Generate first calls via thinning method for non-homog poisson process
        for j=1:nCallTypes
            done = 0;
            lambdaMax = max(arrivalRates(j,:));
            nextTime = 0;
            while done == 0
                nextTime = nextTime + exprnd(60/lambdaMax);
                u = rand;
                if(u <=  arrivalRates(j,min(36,1+floor(nextTime/15)))/lambdaMax) || (nextTime > maxTime)
                    done = 1;
                end
            end
            Events(:,j) = [1; nextTime; j; 0]; % Generate patience times later
        end
        
        %Generate patience times
        RandStream.setGlobalStream(PatientStream);
        for j=1:nCallTypes
            Events(4,j) = Events(2,j)+exprnd(patMean);
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
                    Events = [Events [ 3; shifts(j,1); count; 1]]; % record all agent starts into events
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
            if (prevTime >= 540 && isempty(Events))
                done = 1;
            end
            prevTime = time; % last executed event's time
            
            if isempty(time) == 0
                curShift = 1+floor(time/15);
            end
            
            %Only accept calls until 540, all other events continue to occur
            %after 540 since calls already in the system must be serviced.
            if( Events(1,nextEvent) == 1 )
                if ( time <= maxTime )
                    % CALL ARRIVES
                    callType = Events(3,nextEvent);
                    nCalls(callType, curShift) = nCalls(callType, curShift) + 1; % number of calls by type and shift
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
%                                     pickedCalls(ceil(time), j) = callType; % record call
                                    sTime = exprnd(meanST(callType, nextGroup)); % generate service time
                                    Agents(j,4) = time + sTime; % update availability
                                    Agents(j,3) = 0;
                                    nAvailable(nextGroup) = nAvailable(nextGroup) - 1;
                                    if lastAgent(nextGroup, Agents(j,2)) == -1 && nAvailable(nextGroup) == 0 % when no agent is labeled and this is the last
                                        lastAgent(nextGroup, Agents(j,2)) = j; % label this agent
                                        lastCalls = [lastCalls, [Events(2:3, nextEvent); nextGroup; Agents(j,2); time]]; % record this call
                                    elseif j == lastAgent(Agents(j,1),Agents(j,2)) % when this is labeled as the last agent
                                        lastCalls = [lastCalls, [Events(2:3, nextEvent); nextGroup; Agents(j,2); time]]; % record this call
                                    end
                                    Events = [Events [2; time+sTime; j; 0]]; % schedule call exit
                                    nextGroup = 0;
                                    callTaken = 1;
                                    %Call answered immediately: success
                                    nLessThan20(callType, curShift) = nLessThan20(callType, curShift) + 1;
%                                     totalPicked = totalPicked + 1;
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
                    lambdaMax = max(arrivalRates(callType,:));
                    arrTime = time;
                    while success == 0
                        arrTime = arrTime + exprnd(60/lambdaMax);
                        u = rand;
                        if( u <=  arrivalRates(callType,min(36,1+floor(arrTime/15)))/lambdaMax )
                            success = 1;
                        end
                    end
                    %Generate Patience time
                    RandStream.setGlobalStream(PatientStream);
                    Events = [Events [1; arrTime; callType; arrTime + exprnd(patMean)]];
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
                    while j <= length(callQueue(1,:)) && callTaken == 0
                        if callQueue(3,j) < time % call expiration time
                            lateCalls = [lateCalls, [callQueue(:,j); Inf; NaN]];
                            callQueue(:,j) = []; %Customer expired
                        else
                            if (R(Agents(a,1),callQueue(2,j)) ~= 1000) % take first feasible call in queue
                                %Agent can take this call and abandonment time has
                                %not expired.
                                %Generate Service Time
                                RandStream.setGlobalStream(ServiceStream);
                                sTime = exprnd(meanST(callQueue(2,j),Agents(a,1)));
                                Agents(a,4) = time + sTime;
                                Agents(a,3) = 0;
                                Events = [Events [2; Agents(a,4); a; 0]]; % schedule call exit
                                
                                %Count call as successful if answered withing
                                %SL time constraint.
                                if( time - callQueue(1,j) <= SLtime/60 )
                                    nLessThan20(callQueue(2,j), 1+floor(callQueue(1,j)/15)) = nLessThan20(callQueue(2,j), 1+floor(callQueue(1,j)/15)) + 1;
                                else
                                    % call cannot be answered in time
                                    lateCalls = [lateCalls, [callQueue(:,j); time; sTime]];
                                end
                                
                                callQueue(:,j) = [];
                                callTaken = 1;
                            else
                                j = j+1;
                            end
                        end
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
                agentNumbers = 0;
                
                while Events(2,nextEvent) == time && Events(1,nextEvent) == 3 % agent arrival
                    a = Events(3,nextEvent); % agent group number
                    Agents(a,3) = 1; % set availability = 1
                    Agents(a,4) = Events(2,nextEvent); % next avail = call exit time
                    nAvailable(Agents(a,1)) = 1 + nAvailable(Agents(a,1)); % record number of available agent in that group
                    %Generate "next break" event
                    nextBreak = Events(4,nextEvent); % "starting subcode" - which is the next break for this agent start?
                    agentNumbers = [agentNumbers a]; % record the available agents' group numbers
                    if ( shifts(Agents(a,2), nextBreak+1) < time ) % time of next break is already happening
                        Events = [Events [4; time; a; nextBreak ]]; % make it happen
                    else
                        Events = [Events [4; shifts(Agents(a,2), nextBreak+1); a; nextBreak ]]; % schedule the next break by the time in shifts
                    end
                    Events(:,nextEvent) = []; % clear this executed event
                    [~, nextEvent] = min(Events(2,:)); % find next soonest event in FEL
                end
                
                % Now assign calls to as many of them as possible by FIFO
                j=2; % because first element in callQueue is the all-zero initialization
                while j <= length(callQueue(1,:)) && length(agentNumbers) > 1 % if there are more than j calls in queue and available agents
                    if (callQueue(3,j) > time) % not expired yet
                        callType = callQueue(2,j);
                        minimum = 1000;
                        minAgent = 1000;
                        for k=2:length(agentNumbers)
                            if R(Agents(agentNumbers(k),1),callType) < minimum % R(group, type) is feasible
                                minimum = R(Agents(agentNumbers(k),1),callType); % min preference number for each call in the queue
                                minAgent = k; % the agent of the min preference number (most preferred)
                            end
                        end
                        
                        if minimum ~= 1000 % any feasible agent in agentNumbers
                            %There was an agent newly available to take the
                            %call.
                            %Generating service time of agent
                            RandStream.setGlobalStream(ServiceStream);
                            sTime = exprnd(meanST(callType,Agents(agentNumbers(k),1)));
                            %Count call as successful if answered within
                            %SL time constraint.
                            if( time - callQueue(1,j) <= SLtime/60 ) % time answered within 20s
                                nLessThan20(callQueue(2,j), 1+floor(callQueue(1,j)/15)) = nLessThan20(callQueue(2,j), 1+floor(callQueue(1,j)/15)) + 1;
                            else % nLessThan20 by call type, index of 15 minute intervals
                                lateCalls = [lateCalls, [callQueue(:,j); time; sTime]]; % call is late
                            end
                            a = agentNumbers(minAgent); % pick the agent who is most preferred by the call
%                             pickedCalls(ceil(callQueue(1,j)), a) = callType; % record call
                            Agents(a,3) = 0;
                            Agents(a,4) = time + sTime; % next available time
                            nAvailable(Agents(a,1)) = nAvailable(Agents(a,1)) - 1;
                            if lastAgent(Agents(a,1),Agents(a,2)) == -1 && nAvailable(Agents(a,1)) == 0 % when no agent is labeled and this is the last one
                                lastAgent(Agents(a,1),Agents(a,2)) = a; % label this agent
                                lastCalls = [lastCalls, [callQueue(1:2,j); Agents(a,1); Agents(a,2); time]]; % record this call
                            elseif a == lastAgent(Agents(a,1),Agents(a,2)) % when this is labeled as the last agent
                                lastCalls = [lastCalls, [callQueue(1:2,j); Agents(a,1); Agents(a,2); time]]; % record this call
                            end
                            callQueue(:,j) = []; % delete the scheduled call from the queue
                            Events = [Events [2; time+sTime; a; 0]]; % schedule call exit
                            agentNumbers(minAgent) = []; % delete this agent from available ones
%                             totalPicked = totalPicked + 1;
                        else
                            j = j+1;
                        end
                    else
                        lateCalls = [lateCalls, [callQueue(:,j); Inf; NaN]];
                        callQueue(:,j) = []; % Customer lost                       
                    end
                end
                %In case a call was taken and is expected to leave before next
                %event (previously found) occurs.
                [~, nextEvent] = min(Events(2,:));
                time = Events(2,nextEvent);
                
            elseif (Events(1,nextEvent) == 4)
                % AGENT GOES ON BREAK. Generate next arrival event unless
                % completion of shift.
                
                a = Events(3, nextEvent);
                breakType = Events(4,nextEvent); % index of break that agent is scheduled to have
                breakLength = [15 30 15];
                if breakType ~= 4 % not ending shift
                    if Agents(a,4) <= time % available at this time
                        returnTime = time + breakLength(breakType);
                        nAvailable(Agents(a,1)) = nAvailable(Agents(a,1)) - 1;
                    else
                        %Already unavailable: assume the on-hand call
                        %finishes
                        returnTime = Agents(a,4) + breakLength(breakType);
                    end
                    Agents(a,4) = returnTime;
                    Agents(a,3) = 0;
                    Events = [Events [3; returnTime; a; breakType+1]]; % schedule agent start after break
                else
                    %Agent completed shift, unless he works in the last time
                    %period, where he must stay until all calls are completed.
                    if shifts(Agents(a,2),5) ~= 540 % does not end at end-of-day
                        if(Agents(a,4) < time) % next available time is now
                            nAvailable(Agents(a,1)) = nAvailable(Agents(a,1)) - 1;
                            Agents(a,4) = shifts(Agents(a,2),5);
                            Agents(a,3) = 0;
                        else
                            Agents(a,4) = Agents(a,4) + 1E-10; % to distinguish from those who were available at 540
                            % Note that Agents(:,4) is used to to see what time
                            % they leave at the end of the day. It is only
                            % necessary to increment it by a little just so
                            % that time ~= Agents(a,4) during his next
                            % completion event.
                        end
                    end
                end
                
                Events(:,nextEvent) = [];
                [~, nextEvent] = min(Events(2,:));
                time = Events(2,nextEvent);
            end
        end
        %Loop through agents to update time at which they left (for those who
        %worked until the end of the day).
        for a=1:length(Agents(:,1))
            if ( shifts(Agents(a,2),5) == 540 )
                Agents(a,4) = max(Agents(a,4),540);
                Agents(a,3) = 0;
            end
        end
%         Output(:,:,i) = nLessThan20./nCalls;
%         overallSL(:,:,i) = sum(sum(nLessThan20))/sum(sum(nLessThan20));
        %Set arrival stream for next day
        RandStream.setGlobalStream(ArrivalStream);
        
        SLAll(i,1) = sum(sum(nLessThan20)) / sum(sum(nCalls));
        % Calculate gradient table
        %[fg, bg] = GradientTablePickedCalls(x, beta, nCallTypes, nCalls, meanST, lateCalls, lastCalls, pickedCalls, R, shifts);
        RandStream.setGlobalStream(ServiceStream)
        [fg, bg] = GradientTableOld( x, beta, nCallTypes, nCalls, meanST, lateCalls, lastCalls, R, shifts );
        forwardGradient = forwardGradient + fg;
        backwardGradient = backwardGradient + bg;
    end
    %Return old Random Number Stream
    RandStream.setGlobalStream(OldStream);
    
    % meanSL = mean(Output,3); % mean SL by call type and 15-min period for each day   
    % cost = sum(sum(CostPerDay( x, R, shifts )));
    cost = [1.0, 1.0, 1.1, 1.1, 1.1, 1.2] * x;
    obj = beta * (SLAll - serviceLevelMin) - cost;
    f = mean(obj);
    stdev = std(obj);
    SL = mean(SLAll);
    %constraint = serviceLevelMin - meanSL; % if this has positive components, soln not feasible
    forwardGradient = forwardGradient / runlength;
    backwardGradient = backwardGradient / runlength;
end
