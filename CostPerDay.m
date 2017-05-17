function [ CostPerDay ] = CostPerDay( x, R, shifts )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
varsigma = 0.2; % cost per agent's skill, used in the cost function
CostPerDay = x;
for group = 1:length(x(:,1))
    for shift = 1:length(x(1,:))
        if x(group, shift) > 0
            nu = sum(R(group,:) ~= 1000); % number of callTypes that can be handled
            len = (shifts(shift,5) - shifts(shift,1))/60;
            CostPerDay(group, shift) = x(group, shift) * (1+(nu-1) * varsigma)*(len/30);
        end
    end
end
end

