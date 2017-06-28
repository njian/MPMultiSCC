function [ x ] = evenly_spread( n, row, col )
% Evenly distributed n over row and col

x = ones(row, col) * floor(n/(row*col));   
for k = 1:(n-sum(sum(x)))
    i = randi(row);
    j = randi(col);
    x(i, j) = x(i, j) + 1;
end

end

