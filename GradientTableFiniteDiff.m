function [ forwardGradient, backwardGradient ] = GradientTableFiniteDiff( x, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts )
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
f = MultiSkillPickedCalls(x, beta, ceil(runlength), seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);

forwardGradient = zeros(nGroups, nShifts); % forward gradient, change in obj with one more agent
backwardGradient = zeros(nGroups, nShifts); % backward gradient, change in obj with one less agent
% gradient = change in SL * beta - change in cost
for i = 1:nGroups
    for q = 1:nShifts
        x_trial = x;
        if x(i,q) >= 1
            x_trial(i,q) = x(i,q) - 1;
            backwardGradient(i,q) = MultiSkillPickedCalls(x_trial, beta, 1, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts) - f;
        else
            backwardGradient(i,q) = -Inf;
        end
        
        x_trial = x;
        x_trial(i,q) = x(i,q) + 1;
        forwardGradient(i,q) = MultiSkillPickedCalls(x_trial, beta, 1, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts) - f;
    end
end

end

