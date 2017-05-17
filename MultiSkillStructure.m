function [minmax d m VarNature VarBds FnGradAvail NumConstraintGradAvail, StartingSol, budget, ObjBd, OptimalSol] = MultiSkillStructure(NumStartingSol, seed)
%function [minmax d m VarNature VarBds FnGradAvail NumConstraintGradAvail, StartingSol] = MultiSkillStructure(NumStartingSol, seed);
% Inputs:
%	a) NumStartingSol: Number of starting solutions required. Integer, >= 0
%	b) seed: Seed for generating random starting solutions. Integer, >= 1
% Return structural information on optimization problem
%     a) minmax: -1 to minimize objective , +1 to maximize objective
%     b) d: positive integer giving the dimension d of the domain
%     c) m: nonnegative integer giving the number of constraints. All
%        constraints must be inequality constraints of the form LHS >= 0.
%        If problem is unconstrained (beyond variable bounds) then should be 0.
%     d) VarNature: a d-dimensional column vector indicating the nature of
%        each variable - real (0), integer (1), or categorical (2).
%     e) VarBds: A d-by-2 matrix, the ith row of which gives lower and
%        upper bounds on the ith variable, which can be -inf, +inf or any
%        real number for real or integer variables. Categorical variables
%        are assumed to take integer values including the lower and upper
%        bound endpoints. Thus, for 3 categories, the lower and upper
%        bounds could be 1,3 or 2, 4, etc.
%     f) FnGradAvail: Equals 1 if gradient of function values are
%        available, and 0 otherwise.
%     g) NumConstraintGradAvail: Gives the number of constraints for which
%        gradients of the LHS values are available. If positive, then those
%        constraints come first in the vector of constraints.
%     h) StartingSol: One starting solution in each row, or NaN if NumStartingSol=0.
%        Solutions generated as per problem writeup
%     i) budget: Column vector of suggested budgets, or NaN if none suggested
%     j) ObjBd is a bound (upper bound for maximization problems, lower
%        bound for minimization problems) on the optimal solution value, or
%        NaN if no such bound is known.
%     k) OptimalSol is a d dimensional column vector giving an optimal
%        solution if known, and it equals NaN if no optimal solution is known.

%   ***************************************
%   *** Written by Danielle Lertola to  ***
%   *** use standard calling and random ***
%   *** number streams  6/8/2012        ***
%   ***      Edited by Bryan Chong      ***
%   ***  bhc34@cornell.edu Nov 3rd 2014 ***
%   ***************************************

% Number of groups and call types is user selectable, but in this version
% we use 2 for both.


numShifts=74;
numPeriods=36;

% %%Small Problem%%
numGroups=2;
nCallTypes=2;

%%%%%%%%  Large Problem  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numGroups=35;
% nCallTypes=20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


minmax = -1; % minimize total cost of hiring agents (+1 for maximize)
d = numGroups*numShifts; % # agents working in group i in i=[1,n] during shift q in Q=[1,74].
m = nCallTypes*numPeriods; % one service level requirement for each call type must be satisfied
VarNature = ones(d, 1); % integer variables
VarBds = ones(d, 1) * [0, inf]; % number of agents 
FnGradAvail = 0; % No derivatives 
NumConstraintGradAvail = 0; % No constraint gradients
budget = [100; 1000];
ObjBd=NaN;
OptimalSol=NaN;

if (NumStartingSol < 0) || (NumStartingSol ~= round(NumStartingSol)) || (seed <= 0) || (round(seed) ~= seed),
    fprintf('NumStartingSol should be integer >= 0, seed must be a positive integer\n');
    StartingSol = NaN;
    
    %%Small Problem%%
else
    if (NumStartingSol == 0),
        StartingSol = NaN;
    elseif (NumStartingSol == 1),
        StartingSol=zeros(numGroups,numShifts);
        StartingSol(1,57)=20;
        StartingSol(1,68)=20;
        StartingSol(2,57)=10;
        StartingSol(2,68)=10;

    else
        StartingSol=zeros(numGroups,numShifts,NumStartingSol);
        
        OurStream = RandStream.create('mlfg6331_64'); % Use a different generator from simulation code to avoid stream clashes
        OurStream.Substream = seed;
        OldStream = RandStream.setGlobalStream(OurStream);
                
        edges=[57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74];
        rand1=ceil(18*rand(1,40,NumStartingSol))+56;
        type1=histc(rand1,edges);
        rand2=ceil(18*rand(1,20,NumStartingSol))+56;
        type2=histc(rand2,edges);
            
        StartingSol(:,57:74,:)=[type1;type2];
        
        RandStream.setGlobalStream(OldStream); % Restore previous stream
    end %if NumStartingSol
    
%%%%%%%  Large Problem  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% else
%     if (NumStartingSol == 0),
%         StartingSol = NaN;
%     elseif (NumStartingSol == 1),
%         StartingSol=zeros(numGroups,numShifts);
%         StartingSol(:,57)=3*ones(numGroups,1);
%         StartingSol(:,68)=3*ones(numGroups,1);
%     else
%       StartingSol=zeros(numGroups,numShifts,NumStartingSol);
%         
%       OurStream = RandStream.create('mlfg6331_64'); % Use a different generator from simulation code to avoid stream clashes
%       OurStream.Substream = seed;
%       OldStream = RandStream.setGlobalStream(OurStream);
%                 
%       edges=[57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74];
%       AgentDist=ceil(18*rand(numGroups,6,NumStartingSol))+56;
%       type=histc(AgentDist,edges,2);    
%       StartingSol(:,57:74,:)=type;
%         
%       RandStream.setGlobalStream(OldStream); % Restore previous stream  
%     end %if NumStartingSol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %if inputs ok