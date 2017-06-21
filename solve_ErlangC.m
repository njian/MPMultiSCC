function [ c ] = solve_ErlangC( arrivalRate, meanST, SLtime )

excess = @(c) ErlangC_wait_time(arrivalRate, meanST, c) - SLtime;
c = floor(arrivalRate * meanST);
eTime = excess(c);
while eTime > 0 || isnan(eTime)
   c = c + 1;
   eTime = excess(c);
end

end

