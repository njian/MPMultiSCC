function [ avgWait ] = ErlangC_wait_time( arrivalRate, meanST, c )
% Return the average waiting time of a M/M/c/inf queue with arrivalRate and
% meanST and c servers. arrivalRate and meanST should be in the same time
% unit.

rho = arrivalRate * meanST / c;
if rho >= 1
    avgWait = Inf;
else
    total = 0;
    for m = 0:c-1
        total = total + (c*rho)^m / factorial(m);
    end
    P0 = 1/total;
    avgWait = ((c*rho)^c * rho / factorial(c) / (1-rho)^2 * P0) / arrivalRate;
end

end

