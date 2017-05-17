1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
55
56
57
58
59
60
61
62
63
64
65
66
67
68
69
70
71
72
73
74
75
76
77
78
79
80
81
82
83
84
85
86
87
88
89
90
91
92
93
94
95
96
97
98
99
100
101
102
103
104
105
106
107
108
109
110
111
112
113
114
115
116
117
118
119
120
121
122
123
124
125
function [f_beta, x_ast, SL_beta, sd_beta, n_sim] = localSearch_x(n, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts)
% Given n and beta, local search for x.
nShifts = size(shifts, 1);
 
% Search Parameters
nChangeMax = 2^floor(log2(nAgentGroups*nShifts/10)); % maximum number of swaps in x when generating trial solution
r = 5; % local random search among this many largest gradient components
max_fail = 10; % maximum number of consecutive fails in local search for x
 
count_x = 0;
nChange = nChangeMax; 
count_failed_x = 0;
n_sim = 0;
% x_trial = x_trial_in;
 
% Initialize
x_trial = zeros(nAgentGroups, nShifts);    
for k = 1:n
    row = randi(nAgentGroups);
    col = randi(nShifts);
    x_trial(row, col) = x_trial(row, col) + 1;
end
[f_beta, SL_beta, sd_beta, forwardGradient, backwardGradient] = MultiSkillPickedCalls(x_trial, beta, runlength, seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
backwardGradient(x_trial==0) = -Inf;
 
% f_beta = f_beta_x;
x_ast = x_trial;
% SL_beta = SL_beta_x;
% sd_beta = sd_beta_x;
STOP = 0;
% stop when either consecutive fails exceed tolerance, or there's
% no room for improvement
while count_failed_x < max_fail && STOP == 0
    count_x = count_x + 1;
    fprintf('--------------------------------------------------------------------- \n');
    fprintf('Local search iteration (%d) of x. \n', count_x);
    fprintf('n = %d, beta = %.2f, nChange = %d. \n', sum(sum(x_trial)), beta, nChange);
 
    if count_failed_x == 0
        % Sort the gradients if new gradients are obtained
        [SortedForwardSave, IndForwardSave] = sort(forwardGradient(:),'descend');
        [SortedBackwardSave, IndBackwardSave] = sort(backwardGradient(:),'descend');
    end
    SortedForward = SortedForwardSave;
    IndForward = IndForwardSave;
    SortedBackward = SortedBackwardSave;
    IndBackward = IndBackwardSave;
    % the positive components in forwardGradient and non-negative
    % components in backwardGradient indicate room for improvement
    npForward = sum(SortedForward>=0);
    npBackward = sum(SortedBackward>=0);
    if npForward <=0 && npBackward <= 0
        STOP = 1; % local maximum
    end
    fprintf('Number of positive forward and backward gradients: %d and %d \n', npForward, npBackward);
     
    % search radius = min(r, number of positive components)
    rB = min(r, sum(SortedBackward > -SortedForward(1)));
    rF = min(r, sum(SortedForward > -SortedBackward(1)));
    q = 0;
    while countChange < nChange && rF > 0 && rB > 0 && STOP == 0 %(npForward > 0 || npBackward > 0)
        q = q+1;
        if q > 100
            rF
            rB
        end
        % Generate random group and shift to add an agent
        randIndexF = randi(rF);
        [pRow, pCol] = ind2sub([nAgentGroups, nShifts], IndForward(randIndexF));
        % Generate random group and shift to remove an agent
        randIndexB = randi(rB);
        [mRow, mCol] = ind2sub([nAgentGroups, nShifts], IndBackward(randIndexB));
 
        x_trial(pRow, pCol) = x_trial(pRow, pCol) + 1;            
        IndForward(randIndexF) = []; % remove the chosen index
        SortedForward(randIndexF) = [];
 
        x_trial(mRow, mCol) = x_trial(mRow, mCol) - 1; 
        IndBackward(randIndexB) = [];
        SortedBackward(randIndexB) = [];
 
        countChange = countChange + 1;
         
        % search radius
        rB = min(r, sum(SortedBackward > -SortedForward(1)));
        rF = min(r, sum(SortedForward > -SortedBackward(1)));
    end
     
    if STOP == 0 %~(STOP ~= 0 && count_x == 1)
        [f_beta_x, SL_beta_x, sd_x, forwardGradient_x, backwardGradient_x] = MultiSkillPickedCalls(x_trial, beta, ceil(runlength), seed, serviceLevelMin, nCallTypes, nAgentGroups, arrivalRates, R, Route, shifts);
        fprintf('Trial solution obj = %.2f, SL = %.2f. q = %d. \n', f_beta_x, SL_beta_x, q);        
        n_sim = n_sim + runlength;
 
        if f_beta_x > f_beta
            fprintf('x step: IMPROVE f_beta from %.2f to %.2f +/- %.2f. \n', f_beta, f_beta_x, sd_x);
%             fprintf('x step: n = %d. \n', sum(sum(x_trial)));
            count_failed_x = 0;
 
            % update last best solution for this beta
            x_ast = x_trial;
            f_beta = f_beta_x;
            sd_beta = sd_x;
            SL_beta = SL_beta_x;
            forwardGradient = forwardGradient_x;
            backwardGradient = backwardGradient_x;
            backwardGradient(x_ast==0) = -Inf; % remove impossible reductions in x
        else
            x_trial = x_ast;
            count_failed_x = count_failed_x + 1;
            fprintf('x step: failed for %d times: this solution %.2f < last best %.2f. \n', count_failed_x, f_beta_x, f_beta);
            % update stepsize
            if count_failed_x == max_fail && nChange > 1
                nChange = ceil(nChange/2);
                count_failed_x = 0;
            end
 
        end
    else
        % First parse found no room for improvement. Need to
        % increase/decrease n.
        fprintf('x step: no room for improvement. Seems like a local optimal. \n');
    end
end % end local search for x
 
end