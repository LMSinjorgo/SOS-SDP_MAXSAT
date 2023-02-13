function [Z,x0] = PsdConeProjection(Z)
% Code accompanying the paper:
% On solving the MAX-SAT using sum of squares
% Lennart Sinjorgo, Renata Sotirov
% Feb 2023
% Tilburg, Netherlands
%
% Projects matrix onto PSD cone via its eigenvalue (EV) decomposition.
% Computing the EV decomposition is done in single precision, for speed up over double
% precision.
n = size(Z,1);

% switch to single precision speed
% rename to eigVec to reserve memory
eigVec = single(Z);
eigVec = eigVec+triu(eigVec,1)'; %should be symmetric

%Now compute eigenvalue decomposition of Z
% reuse memory for eigVec
[eigVec, lambda] = eig(eigVec,'vector');

% switch to double for functionality
lambda = double(lambda);

% minimum eigenvector (used for LOBPCG algorithm warm start)
x0 = eigVec(:,1);

% get index where eigenvalues switch from negative to positive
boolIdx = find(lambda >= 10^-8,1,'first');

if isempty(boolIdx)
    Z = zeros(size(Z));
    return;
end
if boolIdx < n/2
    boolCounter = boolIdx - 1;
    boolIdx = 1:boolCounter;
    % switch to double for functionality
    negEigVec = double(eigVec(:,1:boolCounter));
    Z = triu(Z - negEigVec * sparse(boolIdx,boolIdx,lambda(boolIdx),boolCounter,boolCounter) * negEigVec');
else
    boolCounter = n-boolIdx+1;
    boolIdx = boolIdx:n;
    PosEigVec = double(eigVec(:,boolIdx));
    Z = triu(PosEigVec *sparse(1:boolCounter,1:boolCounter,lambda(boolIdx),boolCounter,boolCounter) * PosEigVec');
end
end