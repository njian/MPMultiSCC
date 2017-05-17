function [ forwardGradient, backwardGradient ] = GradientTablePickedCalls( x, beta, nCallTypes, nCalls, meanST, lateCalls, lastCalls, pickedCalls, R, shifts )
%----------------********************************--------------------------
% GradientTable.m
% Description: Gives the estimated change in the objective value 
%              max_x beta * (SL_x - serviceLevelMin) - cost
%              with regard to one more (forwardGradient) or one less
%              (backwardGradient) number of agents in each group * shift.
%   
% Inputs:
%   lateCalls: calls that are picked up late (>20s) or expired without 
%              being picked up. Each column is call arrival time, 
%              call type, call expiry time, pick up time or Inf (if never
%              picked up), service time or NaN (if never picked up)
%   pickedCalls: callType if in minute i agent j picked up a call of
%                this callType else 0. Agent ordered by group and within 
%                group by shift (group 1 shift 1, group 1 shift 2, etc).
%   Others see MultiSkill.m
% Outputs:
%   forwardGradient (backwardGradient): change in the objective wrt/
%             one more (one less) agent in each component of x.
%----------------********************************--------------------------

[nGroups, nShifts] = size(x);
forwardGradient = NaN * ones(nGroups, nShifts); % forward gradient, change in obj with one more agent
backwardGradient = NaN * ones(nGroups, nShifts); % backward gradient, change in obj with one less agent

decreaseLessThan20 = zeros(nCallTypes, 36);
increaseLessThan20 = zeros(nCallTypes, 36);
% gradient = change in SL * beta - change in cost
for group = 1:nGroups
    for shift = 1:nShifts
        % cost of hire one more agent in group i shift q
        hireCost = CostPerDay( [group, shift], R, shifts );
        
        %% Forward gradient
        % go through lateCalls to assign calls to the extra agent
        nextAvail = shifts(shift, 1);
        nextBreak = 1; % 1 and 3: coffee break (15min), 2: lunch break (30min)
        breakLength = [15, 30, 15, 100];
        for j = 1:size(lateCalls,2)
            % check calls in lateCalls that are is not expired or picked up
            % late yet
            if nextBreak < 5 && lateCalls(3,j) > nextAvail && lateCalls(4,j) > nextAvail && R(group, lateCalls(2,j)) < 1000
                % available after mean service time
                if ~isnan(lateCalls(5,j))
                    completeTime = lateCalls(1,j) + lateCalls(5,j);
                else
                    completeTime = lateCalls(1,j) + exprnd(meanST);
                end
                
                if completeTime < shifts(shift, nextBreak)
                    nextAvail = completeTime;
                else
                    nextAvail = completeTime + breakLength(nextBreak); % finishes call and go on break
                    nextBreak = nextBreak + 1;
                end
                % calls picked up in less than 20s is increased
                increaseLessThan20(lateCalls(2,j), 1+floor(lateCalls(1,j)/15)) = increaseLessThan20(lateCalls(2,j), 1+floor(lateCalls(1,j)/15)) + 1;
            end
        end
        changeForwardSL = sum(sum(increaseLessThan20))/sum(sum(nCalls));
        forwardGradient(group, shift) = beta * changeForwardSL - hireCost; 
        
        %% Backward gradient
        % unpick the calls by the "last" agent in each group nPerGroup = sum(schedule,2); 
        lastAgent = sum(sum(x(1:group-1,:))) + sum(x(group,1:shift));
        if x(group, shift) > 0 && lastAgent > 0     
            for k = 1:540 % every period where last agent picked up call
                if pickedCalls(k,lastAgent) ~= 0
                    decreaseLessThan20(pickedCalls(k,lastAgent), 1+floor((k-1)/15)) = decreaseLessThan20(pickedCalls(k,lastAgent), 1+floor((k-1)/15)) - 1;
                end
            end
        end
            changeBackwardSL = sum(sum(decreaseLessThan20))/sum(sum(nCalls));
            backwardGradient(group, shift) = beta * changeBackwardSL + hireCost;
    end
end

end


