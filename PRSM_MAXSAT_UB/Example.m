% Code accompanying the paper:
% On solving the MAX-SAT using sum of squares
% Lennart Sinjorgo, Renata Sotirov
% Feb 2023
% Tilburg, Netherlands

% Compute a SOS-SDP based upper bound the MAX-SAT solution.
clc
clear

% Compute a bound for the included MAX-SAT instance.
% More instances are available at
% http://maxsat.ia.udl.cat/benchmarks/
satName = 's3v70c700-1.cnf';

% directory of the CNF file
satDir = pwd;

% convert the CNF file into a matrix
A = cnfConverter(satDir,satName);

% choose monomial basis x. If left empty, the parsing algorithm will 
% choose a basis according to SOS_p^Q procedure (see paper)
x = [];

% set the Q parameter referring to SOS_p^Q
Q = 50;

% Parse the instance.
PRSMinput = SOS_p_Parse(A,x,Q);

% Specify lower bound (LB), as taken from MSE-2016, 
% see  http://maxsat.ia.udl.cat/detailed/complete-ms-random-table.html
% Note: 700-21 = 679
LB = 679; 

% specify the maximum number of PRSM iterations, or time spent (or both).
maxIter = 600;
maxTime = [];

% specify convergence parameter epsilon (see PRSM.m for precise definition)
epsilon = 0.0001;

% set printing boolean. True for diagnostics, false for no intermediate
% output.
printing = true;

% compute the upper bound (UB). See the function 'PRSM' for more details on parameters.
UB = PRSM(PRSMinput, LB, maxIter, maxTime, epsilon, printing);
