function [ forwardGradient, backwardGradient, lateCalls, lastCalls ] = GradientTablePots( beta, nGroups, totalCalls, meanST, lateCalls, lastCalls, R, costByGroup, updateCoord )
%----------------********************************--------------------------
% GradientTable.m
% Description: Gives the estimated change in the objective value 
%              max_x beta * (SL_x - serviceLevelMin) - cost
%              with regard to one more (gradient_f) or one less
%              (gradient_b) number of agents in each group.
% If updateCoord is empty, do the above. Otherwise return the update in
% gradients and late and last calls list after updates in updateCoord.
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
if isempty(updateCoord)
    nLate = size(lateCalls,2);
    nLast = size(lastCalls,2);
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
    for k = 1:nLast
        backwardGradient(lastCalls(6,k),1) = backwardGradient(lastCalls(6,k),1) - 1/totalCalls; % this group and shift
    %     lastCalls(:,k) = [];
    end
    backwardGradient = beta * backwardGradient + costByGroup;
else
    % coordinates of added agent
    i_p = updateCoord(1); % group
    s_p = updateCoord(2); % shift
    % coordinates of removed agent
    i_m = updateCoord(3); % group
    s_m = updateCoord(4); % shift
    
    % change in backwardGradient
    backwardGradient = - sum(lateCalls(6,:) == i_m)/totalCalls;  
    % For each last call of the removed agent, move it to the lateCalls,
    % because it can no longer be picked up and needed to be queued.
    nLast = size(lastCalls,2);
    for k = 1:nLast
        if lastCalls(6,k) == i_m;
           lateCalls = [lateCalls, lastCalls(1:5,k)];
           lastCalls(:,k) = []; 
        end
    end
    
    nLate = size(lateCalls,2);
    increaseLessThan20 = 0;
    nextAvail = 0;
    j = 1;
    % for every late call that can be picked up by group i
    while j <= nLate && R(i_p, lateCalls(2,j)) < 1000
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
            lateCalls(:,j) = [];
        end
        j = j + 1;
    end
    changeForwardSL = increaseLessThan20/totalCalls;
    forwardGradient = beta * changeForwardSL - costByGroup(i_p,1);    
end

end

