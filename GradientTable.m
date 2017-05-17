function [ forwardGradient, backwardGradient ] = GradientTable( x, beta, nCallTypes, nCalls, meanST, lateCalls, lastCalls, R, shifts )
%----------------********************************--------------------------
% GradientTable.m
% Description: Gives the estimated change in the objective value 
%              max_x beta * (SL_x - serviceLevelMin) - cost
%              with regard to one more (gradient_f) or one less
%              (gradient_b) number of agents in each group * shift.
%   
% Inputs:
%   lateCalls: calls that are picked up late (>20s) or expired without 
%              being picked up. Each column is call arrival time, 
%              call type, call expiry time, pick up time or Inf (if never
%              picked up), service time or NaN (if never picked up)
%   lastCalls: calls that are picked up by an agent who is the only 
%              available one in the group. Each column is the same as
%              callQueue, plus the agent group number.
% Outputs:
%   gradient: change in the overall service level with regard to change in
%             each component in x
%----------------********************************--------------------------

[nGroups, nShifts] = size(x);
forwardGradient = zeros(nGroups, nShifts); % forward gradient, change in obj with one more agent
backwardGradient = zeros(nGroups, nShifts); % backward gradient, change in obj with one less agent

decreaseLessThan20 = zeros(nCallTypes, 36);
increaseLessThan20 = zeros(nCallTypes, 36);
% gradient = change in SL * beta - change in cost
for i = 1:nGroups
    for q = 1:nShifts
        % cost of hire one more agent in group i shift q
        hireCost = CostPerDay( [i, q], R, shifts );
        
        %% Forward gradient
        % go through lateCalls to assign calls to the extra agent
        nextAvail = shifts(q, 1);
        nextBreak = 1; % 1 and 3: coffee break (15min), 2: lunch break (30min)
        breakLength = [15, 30, 15, 100];
        for j = 1:size(lateCalls,2)
            % check calls in lateCalls that are is not expired or picked up
            % late yet
            if nextBreak < 5 && lateCalls(3,j) > nextAvail && lateCalls(4,j) > nextAvail && R(i, lateCalls(2,j)) < 1000
                % available after mean service time
                if ~isnan(lateCalls(5,j))
                    completeTime = lateCalls(1,j) + lateCalls(5,j);
                else
                    completeTime = lateCalls(1,j) + exprnd(meanST);
                end
                
                if completeTime < shifts(q, nextBreak)
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
        forwardGradient(i, q) = beta * changeForwardSL - hireCost; 
        
        %% Backward gradient
        % unpick the calls by the "last" agent in each group
        for k = 1:size(lastCalls, 2)
            decreaseLessThan20(lastCalls(4,k), 1+floor(lastCalls(1,k)/15)) = decreaseLessThan20(lastCalls(4,k), 1+floor(lastCalls(1,k)/15)) - 1;
        end
        changeBackwardSL = sum(sum(decreaseLessThan20))/sum(sum(nCalls));  
        backwardGradient(i, q) = beta * changeBackwardSL + hireCost;
    end
end

end

