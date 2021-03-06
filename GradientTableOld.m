function [ forwardGradient, backwardGradient ] = GradientTableOld( x, beta, totalCalls, meanST, lateCalls, lastCalls, R, shifts )
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
%              available one in the group. Each column is call arrival, 
%              call type, agent group, shift, and pick up time.
% Outputs:
%   gradient: change in the overall service level with regard to change in
%             each component in x
%----------------********************************--------------------------

[nGroups, nShifts] = size(x);
breakLength = [15, 30, 15, 1000];
hireCostAll = CostPerDay( ones(nGroups, nShifts), R, shifts );

forwardGradient = zeros(nGroups, nShifts); % forward gradient, change in obj with one more agent
% gradient = change in SL * beta - change in cost
if ~isempty(lateCalls)
    for i = 1:nGroups
        for q = 1:nShifts
%             decreaseLessThan20 = zeros(nCallTypes, 36);
            increaseLessThan20 = 0;
            
            % cost of hire one more agent in group i shift q
            hireCost = hireCostAll(i, q);

            %% Forward gradient
            % go through lateCalls to assign calls to the extra agent
            nextAvail = shifts(q, 1);
            nextBreak = 1; % 1 and 3: coffee break (15min), 2: lunch break (30min)
            for j = 1:size(lateCalls,2)
                % check calls in lateCalls that are is not expired or picked up
                % late yet
                if nextBreak < 5 && lateCalls(3,j) > nextAvail && lateCalls(4,j) > nextAvail && R(i, lateCalls(2,j)) < 1000
                    % available after mean service time
                    if ~isnan(lateCalls(5,j))
                        completeTime = lateCalls(1,j) + lateCalls(5,j);
                    else
                        completeTime = lateCalls(1,j) + exprnd(meanST(lateCalls(2,j),i));
                    end

                    if completeTime < shifts(q, nextBreak+1)
                        nextAvail = completeTime;
                    else
                        nextAvail = completeTime + breakLength(nextBreak); % finishes call and go on break
                        nextBreak = nextBreak + 1;
                    end
                    % calls picked up in less than 20s is increased
                    increaseLessThan20 = increaseLessThan20 + 1;
                end
            end
            changeForwardSL = increaseLessThan20/totalCalls;
            forwardGradient(i, q) = beta * changeForwardSL - hireCost; 
        end
    end
else
    forwardGradient = - hireCostAll; % no late calls
end

%% Backward gradient
% unpick the calls by the "last" agent in each group
backwardGradient = - beta * lastCalls / totalCalls + hireCostAll;

end

