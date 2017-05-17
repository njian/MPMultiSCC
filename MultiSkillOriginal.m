function [nLessThan20, nCalls] = MultiSkillOriginal(x, runlength, seed, ~)
% function [fn, FnVar, FnGrad, FnGradCov, constraint, ConstraintCov, ConstraintGrad, ConstraintGradCov] = MultiSkill(x, runlength, seed, other);
% x is a matrix (nGroups x nShifts) containing the number of agents in
% group i working shift j
% runlength is the number of days of simulated time to simulate
% seed is the index of the substreams to use (integer >= 1)
% other is not used
% RETURNS: Mean and confidence interval of the proportion of deliveries
% completed in less than "Tau".
% Uses discrete event simulation with only 2 events:
% 1. Order Arrives
% 2. Truck arrives back to warehouse.
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
%   ***************************************

FnVar = NaN;
FnGrad = NaN;
FnGradCov = NaN;
ConstraintCov = NaN;
ConstraintGrad = NaN;
ConstraintGradCov = NaN;



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

%Text file shifts contains all of the shifts being considered in the
%following format:
% col 1. Time shift starts, in minutes.
% col 2. Time of first break (15 min).
% col 1. Time of lunch break (30 min).
% col 2. Time of second break (15 min).
% col 1. Time shift ends
shifts = csvread('Shifts.csv');

nCallTypes =2;
nAgentGroups = 2;
varsigma = 0.2;
arrivalRates = csvread('ArrivalRatesSmall.csv');
serviceLevelMin = 0.8;
meanST = 8;
patMean = 10;
SLtime = 20;
nDays = runlength;


if (sum(sum(x < 0)) > 0|| (sum(size(x) ~= [nCallTypes, 74])>0) || (runlength <= 0) || (seed <= 0) || (round(seed) ~= seed))
    fprintf('x should be %d * 74 matrix with nonnegative components\nrunlength should be positive and real,\nseed should be a positive integer\n', nCallTypes);
    fn = NaN;
    FnVar = NaN;
else % main simulation
    
    % Defines the different agent groups and their preferences. In this case,
    % there are two groups. the first one can only take calls of type 1 while
    % the second one can take both types, yet it prefers type 2 calls. A value
    % of 1000 means that calls cannot be taken. (Inf in problem statement, but
    % 1000 is easier to handle in excel)
    R = [1 1000; 2 1]; % nGroups x nCallTypes
    
    % I x (number of shifts) matrix with the number of agents
    % in each group and schedule being hired.
    schedule = x;
    
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
    Rcopy = R;
    Route = zeros(nCallTypes, nAgentGroups);
    for j=1:nCallTypes
        [C,index] = min(Rcopy(:,j));
        count = 1;
        while(C ~= 1000)
            Route(j,count) = index;
            count = count+1;
            Rcopy(index,j)=1000;
            [C,index] = min(Rcopy(:,j));
        end
    end
    clear Rcopy;
    
    maxTime = 540;
    Output = zeros(nCallTypes, 36, nDays);
    
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
        
        Events = zeros(4,nCallTypes);
        callQueue = zeros(3,1); % (time call entered queue; call type; call expiry time), equivalent to Events(2:4,i) for an arrival event.
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
        Agents = zeros(sum(sum(schedule)), 4);
        count = 1;
        nAvailable = zeros(nAgentGroups,1);
        s = size(schedule);
        for l=1:s(1)         % Agent Group
            for j=1:s(2)       % Agent Shift
                for k=1:schedule(l,j)
                    Agents(count,:)= [l, j, 0, inf];
                    %Type 3: Agent starts: type 1
                    Events = [Events [ 3; shifts(j,1); count; 1]];
                    count = count + 1;
                end
            end
        end
        
        [~, nextEvent] = min(Events(2,:));
        time = Events(2,nextEvent);
        done = 0;
        prevTime = 0;
        while( done == 0 )
            %Check if all events have been dealt with, time >= 540
            if (prevTime >= 540 && isempty(Events))
                done = 1;
            end
            prevTime = time;
            
            if isempty(time) == 0
                curShift = 1+floor(time/15);
            end
            
            %Only accept calls until 540, all other events continue to occur
            %after 540 since calls already in the system must be serviced.
            if( Events(1,nextEvent) == 1 )
                if ( time <= maxTime )
                    % CALL ARRIVES
                    callType = Events(3,nextEvent);
                    nCalls(callType, curShift) = nCalls(callType, curShift) + 1;
                    count = 1;
                    nextGroup = Route(callType,count);
                    callTaken = 0;
                    while nextGroup ~= 0 && count <= nAgentGroups
                        if nAvailable(nextGroup) > 0
                            %Set stream for generating service time of agent
                            RandStream.setGlobalStream(ServiceStream);
                            %find agent, takes the call
                            for j = 1+sum(nPerGroup(1:nextGroup-1)):sum(nPerGroup(1:nextGroup))
                                if( Agents(j,4) <= time && Agents(j,3) == 1 )
                                    sTime = exprnd(meanST);
                                    Agents(j,4) = time + sTime;
                                    Agents(j,3) = 0;
                                    nAvailable(nextGroup) = nAvailable(nextGroup) - 1;
                                    Events = [Events [2; time+sTime; j; 0]];
                                    nextGroup = 0;
                                    callTaken = 1;
                                    %Call answered immediately: success
                                    nLessThan20(callType, curShift) = nLessThan20(callType, curShift) + 1;
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
                
                a = Events(3,nextEvent);
                if Agents(a,4) == time
                    %Agent becomes available: Search for a call in queue and,
                    %if none, set Agent as idle.
                    j=2; callTaken = 0;
                    while j <= length(callQueue(1,:)) && callTaken == 0
                        if callQueue(3,j) < time
                            callQueue(:,j) = []; %Customer expired
                        else
                            if (R(Agents(a,1),callQueue(2,j)) ~= 1000)
                                %Agent can take this call and abandonment time has
                                %not expired.
                                %Generate Service Time
                                RandStream.setGlobalStream(ServiceStream);
                                sTime = exprnd(meanST);
                                Agents(a,4) = time + sTime;
                                Agents(a,3) = 0;
                                Events = [Events [2; time+sTime; a; 0]];
                                
                                %Count call as successful if answered withing
                                %SL time constraint.
                                if( time - callQueue(1,j) <= SLtime/60 )
                                    nLessThan20(callQueue(2,j), 1+floor(callQueue(1,j)/15)) = nLessThan20(callQueue(2,j), 1+floor(callQueue(1,j)/15)) + 1;
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
                        Agents(a,3) = 1;
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
                
                while Events(2,nextEvent) == time && Events(1,nextEvent) == 3
                    a = Events(3,nextEvent);
                    Agents(a,3) = 1;
                    Agents(a,4) = Events(2,nextEvent);
                    nAvailable(Agents(a,1)) = 1 + nAvailable(Agents(a,1));
                    %Generate "next break" event
                    nextBreak = Events(4,nextEvent);
                    agentNumbers = [agentNumbers a];
                    if ( shifts(Agents(a,2), nextBreak+1) < time )
                        Events = [Events [4; time; a; nextBreak ]]; % doesn't really happen with params that make sense
                    else
                        Events = [Events [4; shifts(Agents(a,2), nextBreak+1); a; nextBreak ]];
                    end
                    Events(:,nextEvent) = [];
                    [~, nextEvent] = min(Events(2,:));
                end
                
                % Now assign calls to as many of them as possible
                j=2;
                while j <= length(callQueue(1,:)) && length(agentNumbers) > 1
                    if (callQueue(3,j) > time)
                        callType = callQueue(2,j);
                        minimum = 1000;
                        minAgent = 1000;
                        for k=2:length(agentNumbers)
                            if R(Agents(agentNumbers(k),1),callType) < minimum
                                minimum = R(Agents(agentNumbers(k),1),callType);
                                minAgent = k;
                            end
                        end
                        
                        if minimum ~= 1000
                            %There was an agent newly available to take the
                            %call.
                            %Count call as successful if answered within
                            %SL time constraint.
                            if( time - callQueue(1,j) <= SLtime/60 )
                                nLessThan20(callQueue(2,j), 1+floor(callQueue(1,j)/15)) = nLessThan20(callQueue(2,j), 1+floor(callQueue(1,j)/15)) + 1;
                            end
                            callQueue(:,j) = [];
                            %Generating service time of agent
                            RandStream.setGlobalStream(ServiceStream);
                            sTime = exprnd(meanST);
                            a = agentNumbers(minAgent);
                            Agents(a,3) = 0;
                            Agents(a,4) = time + sTime;
                            nAvailable(Agents(a,1)) = nAvailable(Agents(a,1)) - 1;
                            Events = [Events [2; time+sTime; a; 0]];
                            agentNumbers(minAgent) = [];
                        else
                            j = j+1;
                        end
                    else
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
                breakType = Events(4,nextEvent);
                breakLength = [15 30 15];
                if breakType ~= 4
                    if Agents(a,4) <= time
                        returnTime = time + breakLength(breakType);
                        nAvailable(Agents(a,1)) = nAvailable(Agents(a,1)) - 1;
                    else
                        %Already unavailable
                        returnTime = Agents(a,4) + breakLength(breakType);
                    end
                    Agents(a,4) = returnTime;
                    Agents(a,3) = 0;
                    Events = [Events [3; returnTime; a; breakType+1]];
                else
                    %Agent completed shift, unless he works in the last time
                    %period, where he must stay until all calls are completed.
                    if shifts(Agents(a,2),5) ~= 540
                        if(Agents(a,4) < time)
                            nAvailable(Agents(a,1)) = nAvailable(Agents(a,1)) - 1;
                            Agents(a,4) = shifts(Agents(a,2),5);
                            Agents(a,3) = 0;
                        else
                            Agents(a,4) = Agents(a,4) + 1E-10;
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
        Output(:,:,i) = nLessThan20./nCalls;
        %Set arrival stream for next day
        RandStream.setGlobalStream(ArrivalStream);
    end
    %Return old Random Number Stream
    RandStream.setGlobalStream(OldStream);
    
    meanSL = mean(Output,3);
    
    CostPerDay = 0;
    for i = 1:length(Agents(:,1))
        nu = sum(R(Agents(i,1),:) ~= 1000);
        len = (shifts(Agents(i,2),5) - shifts(Agents(i,2),1))/60;
        CostPerDay = CostPerDay +(1+(nu-1) * varsigma)*(len/30);
    end
    
    fn = CostPerDay;
    constraint = serviceLevelMin - meanSL; % if this has positive components, soln not feasible
    
end
