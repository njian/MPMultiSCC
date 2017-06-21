function [ avgWait ] = ErlangC_wait_time( arrivalRate, meanST, c )
% Return the average waiting time of a M/M/c/inf queue with arrivalRate and
% meanST and c servers

rho = arrivalRate * meanST / c;

total = 0;
for m = 0:c-1
    total = total + (c*rho)^m / factorial(m);
end
total = total + (c*rho)^c / factorial(c) / (1-rho);

avgWait = ((c*rho)^c * rho / factorial(c) / (1-rho)^2 / total) / arrivalRate;

end

