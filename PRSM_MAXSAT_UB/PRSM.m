function [UB, Z,M,S, stoppingReason] = PRSM(PRSMinput, LB, maxIter, maxTime, epsilon, printing)
% Code accompanying the paper:
% On solving the MAX-SAT using sum of squares
% Lennart Sinjorgo, Renata Sotirov
% Feb 2023
% Tilburg, Netherlands
%
% compute an SDP-SOS based upper bound on the MAX-SAT solution
%
% Inputs:
% PRSMinput: struct, obtained from the parsing function.
% LB: integer. Lowerbound on the MAX-SAT, can be left empty.
% maxIter: integer > 0. maximum number of PRSM iterations to run.
% maxTime: number > 0. maximum time to spend on the PRSM iterations.
% epsilon: number > 0 (convergence criteria). If (objective_{k+1} -
% objective_{k} < epsilon) terminate the PRSM. (Do not put this value too low, e.g.
% epsilon = 10^(-8)).
% printing: boolean. if true, prints diagnostics.
% 
% All inputs except PRSMinput can be empty, or not inputted at all.
%
% Outputs:
% UB: valid upper bound to the MAX-SAT. (i.e., it is impossible to satisfy
% more than UB number of clauses)
% Z,M,S: matrix variables of the PRSM.
% stoppingReason: string explaining the reason for stopping the PRSM.
tStart = tic;
if nargin <= 2
    % defaulting to standard inputs
    maxIter = Inf;
    maxTime = Inf;
    epsilon = 0.0001;
    printing = true;
    if nargin == 1
        LB = [];
    end
else
    if nargin ~= 6
        error("incorrect inputs")
    end
    if isempty(maxIter) && isempty(maxTime) && isempty(epsilon)
        error("No valid stopping criteria given");
    end
    if isempty(maxIter)
        maxIter = Inf;
    end
    if isempty(maxTime)
        maxTime = Inf;
    end
    if isempty(printing)
        printing = true;
    end
end



% get the inputs from the PRSM input struct
constant = PRSMinput.constant; % constant term of the polynomial F_phi

coeffs = PRSMinput.coeffs; % coefficients of F_phi
numClauses = PRSMinput.numClauses; % number of clauses in the MAX-SAt instance

groupIndicators = PRSMinput.groupIndicators;
nonUnitIdx = PRSMinput.nonUnitIdx;
A = [groupIndicators,nonUnitIdx];
nonUnitIdx = sortrows(A);

nonUnitIdx = nonUnitIdx(:,2);


groupSizes = PRSMinput.groupSizes;

groupIndicators = repelem(1:uint32(size(groupSizes,1)),groupSizes)';

invGroupSizes = double(groupSizes) .^ (-1);

% size of the basis
s = size(PRSMinput.x,1);

%Get linear diagonal indices
diagIdx = (1:1+s:s^2)';

% fixed stepsize parameters: gamma1, gamma2, beta
gamma1 = 0.5;
gamma2 = (1+sqrt(5))/2*0.9;
beta = s/2.5;

% Initialize all matrices as 0 matrix.
Z = zeros(s,s); M = Z; Q = Z;

objCorrection = numClauses - constant;

% counts the iteration number
t = 1;
objCheckCounter = 1;

k = 20; % check objective value every k iterations


stoppingCriteriaMet = 0;
stoppingReason = "";


% project matrix M onto set M_phi, as a warm start
projM = M;
M = diag(projM(diagIdx));
nonUnitValues = projM(nonUnitIdx);
subtraction = (accumarray(groupIndicators,nonUnitValues) - coeffs) .* invGroupSizes;
subtractionREPMAT = subtraction(groupIndicators);
M(nonUnitIdx) = nonUnitValues - subtractionREPMAT;


if printing
    fprintf("Matrix size: %d \n", s);
    if ~isempty(LB)
        fprintf("LB: %d \n", LB);
    else
        fprintf("No lower bound given \n")
    end
    fprintf("Maximum num. Iterations: %d \n", maxIter);
    fprintf("Maximum time: %.2f \n", maxTime);
    fprintf("Epsilon: %f \n", epsilon);
    fprintf("(gamma1, gamma2, beta) = (%.2f, %.2f, %.2f) \n", gamma1,gamma2,beta);
    fprintf('------------------- \n');
    fprintf('Iter.    Time    UB \n')
end

if isempty(LB)
    LB = -99;
end

% note that matrices Z,M,S are stored as upper triangular, to skip out
% on the Z = Z+triu(Z,1)' operation
UB_old = Inf;
objDiffPerIter = Inf;
while stoppingCriteriaMet == 0

    % project {M + (1/beta)*S} = {M+Q} onto the PSD cone
    [Z,x0] = PsdConeProjection(M+Q);

    % Update S^k to S^{k+0.5} (via Q)
    Q = Q + gamma1  * (M - Z);

    % M subproblem
    % project {Z - (1/beta)(S+I)} = {Z-Q} onto M_phi
    % adjust diagonal later
    projM = Z - Q;

    % adjust diagonal    
    M(diagIdx) = projM(diagIdx)-(1/beta);
    nonUnitValues = projM(nonUnitIdx);
    subtraction = (accumarray(groupIndicators,nonUnitValues) - coeffs) .* invGroupSizes;
    M(nonUnitIdx) = nonUnitValues - subtraction(groupIndicators);
    
    % Update S^k to S^{k+1} (via Q)    
    Q = Q + gamma2 *(M-Z);

    % compute estimate of upper bound
    UB = sum(M(diagIdx));

    if objCheckCounter >= k || t >= maxIter || toc(tStart) > maxTime
        
        if t < maxIter && toc(tStart) < maxTime
            % compute eigenvalue shift with LOBPCG alg. if this is not final iteration
            [~,minEig] = lobpcg(x0,single(M+triu(M,1)'));
            UB = UB + objCorrection - minEig * s;
        else
            % compute eigenvalue shift with MATLAB eig for final iteration
            UB = LB-1;
        end

        if UB <= LB + 1
            % check if eigenvalue computation was accurate
            minEig = eig(M+triu(M,1)', 'vector');
            minEig = minEig(1);
            UB = sum(M(diagIdx)) + objCorrection - minEig(1) * s;
            if floor(UB) <= LB
                if printing
                    secs_clock = toc(tStart);
                    fprintf('%5.0d %8.3f   %4.5f \n', t, secs_clock, UB);
                end
                stoppingReason = "Objective reached LB";
                break;
            end 
        end
        if printing
            secs_clock = toc(tStart);
            fprintf('%5.0d %8.3f   %4.5f \n', t, secs_clock, UB);
        end

        % compute objective decrease per iteration, for stopping criteria
        objDiffPerIter = (UB_old - UB)/k;
        UB_old = UB;
        
        objCheckCounter = 0; % reset objCheckCounter
    end

    % check stopping criteria
    if objDiffPerIter > 0 && objDiffPerIter < epsilon
        stoppingReason = "Objective has converged";
        break;
    end
    if t >= maxIter || toc(tStart) > maxTime
        if t >= maxIter
            stoppingReason = "Maximum number of iterations";
        else
            stoppingReason = "Maximum time spent";
        end
        break;
    end
    t = t+1;
    objCheckCounter = objCheckCounter + 1;
end


if printing
    fprintf(stoppingReason);
    fprintf('\n------------------- \n');
end

% make M feasible by shifting the eigenvalues
M(diagIdx) = M(diagIdx)-minEig;

% retrieve S for output
S = beta*Q;

end