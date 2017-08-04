function [ forwardGradient, backwardGradient ] = GradientTableStaffing( beta, nGroups, totalCalls, meanST, lateCalls, lastCalls, R, costByGroup )
%----------------********************************--------------------------
% GradientTable.m
% Description: Gives the estimated change in the objective value 
%              max_x beta * (SL_x - serviceLevelMin) - cost
%              with regard to one more (gradient_f) or one less
%              (gradient_b) number of agents in each group.
%   
% Inputs:
%   lateCalls: calls that are picked up late (>20s) or expired without 
%              being picked up. Each column is call arrival time, 
%              call type, call expiry time, pick up time or Inf (if never
%              picked up), service time or NaN (if never picked up)
%   lastCalls: number of calls picked up by the "last" agent in the group.
% Outputs:
%   gradient: change in the overall service level with regard to change in
%             each component in x
%----------------********************************--------------------------
nLate = size(lateCalls,2);
forwardGradient = zeros(nGroups, 1); % forward gradient, change in obj with one more agent
backwardGradient = zeros(nGroups, 1); % backward gradient, change in obj with one less agent
% gradient = change in SL * beta - change in cost
if ~isempty(lateCalls)
    for i = 1:nGroups
        increaseLessThan20 = 0;
        %% Forward gradient
        % go through lateCalls to assign calls to the extra agent
        nextAvail = 0;
        j = 1;
        % for every late call that can be picked up by group i
        while j <= nLate && R(i, lateCalls(2,j)) < 1000
            % check calls in lateCalls that are is not expired or picked up
            % late yet (still in queue)
            if lateCalls(3,j) > nextAvail && lateCalls(4,j) > nextAvail
                % available after mean service time
                if ~isnan(lateCalls(5,j))
                    nextAvail = lateCalls(1,j) + lateCalls(5,j);
                else
                    nextAvail = lateCalls(1,j) + exprnd(meanST(lateCalls(2,j),i));
                end
                % calls picked up in less than 20s is increased
                increaseLessThan20 = increaseLessThan20 + 1;
            end
            j = j + 1;
        end
        changeForwardSL = increaseLessThan20/totalCalls;
        forwardGradient(i, 1) = beta * changeForwardSL - costByGroup(i,1); 
    end
else
    forwardGradient = -costByGroup; % no late calls
end

%% Backward gradient
% unpick the calls by the "last" agent in each group
for i = 1:nGroups
    backwardGradient(i,1) = backwardGradient(i,1) - lastCalls(i,1)/totalCalls;
end
backwardGradient = beta * backwardGradient + costByGroup;
end

