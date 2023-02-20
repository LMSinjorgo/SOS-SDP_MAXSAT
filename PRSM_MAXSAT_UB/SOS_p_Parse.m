function PRSMinput = SOS_p_Parse(Q,x,extraQuadNum)
% Code accompanying the paper:
% On solving the MAX-SAT using sum of squares
% Lennart Sinjorgo, Renata Sotirov
% Feb 2023
% Tilburg, Netherlands
%
% Parses the CNF instance in a format suited for the PRSM.
% For use instructions, see Example.m file.

n = size(Q,2)/2;
m = size(Q,1);
A = Q(:,1:n) - Q(:,(n+1):2*n);

%%%%%% compute F_phi
ell = full(sum(abs(A),2));
constant = sum(0.5 .^ ell);
quadCounter = zeros(n,n);
linearCounter = zeros(n,1);
%%% check all the unique length 3 clauses
absA = abs(A);
[i,j] = find(absA);
absA = sparse(i,j,true,m,n);

% find unique rows in absA
[sortA,~] = sortrows(absA);

groupsSortA = sortA(1:m-1,:) ~= sortA(2:m,:);
groupsSortA = any(groupsSortA,2);
groupsSortA = [true; groupsSortA];
uniqueVariableClauses = sortA(groupsSortA,:);

ellUnique = full(sum(uniqueVariableClauses,2));
unique3clauses = uniqueVariableClauses(ellUnique == 3,:);
unique4clauses = uniqueVariableClauses(ellUnique == 4,:);

numUniq3Clauses = size(unique3clauses,1);
numUniq4Clauses = size(unique4clauses,1);
% cubicCounter checks all the possible cubics appearing in F_phi
% the first 3 columns refer to the indices of the monomials
% the last column takes the coefficient
cubicCounter = find(unique3clauses') - repelem(0:n:n*numUniq3Clauses-1,3)';
quarticCounter = find(unique4clauses') - repelem(0:n:n*numUniq4Clauses-1,4)';

% create monomial basis first inbetween
cubicMapping = [1;2;3;1;2;4;1;3;4;2;3;4];
if isempty(x)
    x = zeros(0,2);
    quadCreation = [1;2;1;3;2;3];
    if any(ellUnique == 2)
        unique2clauses = uniqueVariableClauses(ellUnique == 2,:);
        numUniq2Clauses = size(unique2clauses,1);
        unique2clauses = find(unique2clauses') - repelem(0:n:n*numUniq2Clauses-1,2)';
        x = [x;unique2clauses];
    end
    if any(ellUnique == 3)
        x = [x;cubicCounter(repmat(quadCreation,numUniq3Clauses,1)+repelem((0:3:3*(numUniq3Clauses-1))',6))];
    end
    if any(ellUnique == 4)
        quadCreation4 = [1;2;1;3;2;3;1;4;2;4;3;4];
        x = [x; quarticCounter(repmat(quadCreation4,numUniq4Clauses,1)+repelem((0:4:4*(numUniq4Clauses-1))',12))];
    end
    x = reshape(x,2,size(x,1)/2)';

    % add extra quadratic terms (SOS_p^Q basis)
    if ~isempty(extraQuadNum)
        occurence = sum(absA,1);
        [~,idx] = mink(occurence,extraQuadNum);
        idx = sort(idx);
        x = [x; nchoosek(idx,2);];
    end
    
    numRows = size(x,1);
    % fast matlab 'unique' operator call
    [sortA,~] = matlab.internal.math.sortrowsHelper(x);
    grpSortA = sortA(1:numRows-1,:) ~= sortA(2:numRows,:);
    grpSortA = [true; any(grpSortA,2)];
    x = sortA(grpSortA,:);

    numRows = size(x,1);
    x = x(numRows:-1:1,:)';
    x = sparse([2:n+1,repelem(n+2:n+1+numRows,2)],[(1:n)';x(:)],true,n+1+numRows,n);
end


cubicCounter = reshape(cubicCounter,3,numUniq3Clauses)';
cubicCounter = [cubicCounter, zeros(numUniq3Clauses,1)];

quarticCounter = reshape(quarticCounter,4,numUniq4Clauses)';
quarticCounter = [quarticCounter, zeros(numUniq4Clauses,1)];
numQuartics = size(quarticCounter,1);

for j = 1:m
    aj = A(j,:);
    vars = find(aj ~= 0);
    linearCounter(vars) = linearCounter(vars) - aj(vars)'/(2^ell(j));

    if ell(j) == 2
        quadCounter(vars(2),vars(1)) = quadCounter(vars(2),vars(1)) + aj(vars(1))*aj(vars(2)) / 4;
    elseif ell(j) == 3
        % divide by 8 for (1/2^3) = 1/8
        quadCounter(vars(2),vars(1)) = quadCounter(vars(2),vars(1)) + aj(vars(1))*aj(vars(2)) / 8;
        quadCounter(vars(3),vars(1)) = quadCounter(vars(3),vars(1)) + aj(vars(1))*aj(vars(3)) / 8;
        quadCounter(vars(3),vars(2)) = quadCounter(vars(3),vars(2)) + aj(vars(2))*aj(vars(3)) / 8;

        % ismember is slow matlab function
        %[~,matchingIdx] = ismember(vars, cubicCounter(:,1:3) ,'rows');
        matchingIdx = repmat(vars,size(cubicCounter,1),1);
        matchingIdx = matchingIdx == cubicCounter(:,1:3);
        matchingIdx = all(matchingIdx,2);
        matchingIdx = find(matchingIdx,1,'first');

        cubicCounter(matchingIdx,4) = cubicCounter(matchingIdx,4) - aj(vars(1))*aj(vars(2))*aj(vars(3)) / 8;
    elseif ell(j) == 4
        % divide by 16 for (1/2^4) = 1/16
        quadCounter(vars(2),vars(1)) = quadCounter(vars(2),vars(1)) + aj(vars(1))*aj(vars(2)) / 16;
        quadCounter(vars(3),vars(1)) = quadCounter(vars(3),vars(1)) + aj(vars(1))*aj(vars(3)) / 16;
        quadCounter(vars(3),vars(2)) = quadCounter(vars(3),vars(2)) + aj(vars(2))*aj(vars(3)) / 16;
        quadCounter(vars(1),vars(4)) = quadCounter(vars(3),vars(2)) + aj(vars(1))*aj(vars(4)) / 16;
        quadCounter(vars(2),vars(4)) = quadCounter(vars(3),vars(2)) + aj(vars(2))*aj(vars(4)) / 16;
        quadCounter(vars(3),vars(4)) = quadCounter(vars(3),vars(2)) + aj(vars(3))*aj(vars(4)) / 16;

        % match the cubic coeffs
        cubics = reshape(vars(cubicMapping),3,4);
        for k = 1:4
            currentCubic = cubics(:,k)';
            newCoeff = -full(prod(aj(currentCubic))) / 16;

            % ismember is slow matlab function
            %[~,matchingIdx] = ismember(vars, cubicCounter(:,1:3) ,'rows');
            matchingIdx = repmat(currentCubic,size(cubicCounter,1),1);
            matchingIdx = matchingIdx == cubicCounter(:,1:3);
            matchingIdx = all(matchingIdx,2);
            matchingIdx = find(matchingIdx,1,'first');
            if isempty(matchingIdx)
                cubicCounter = [cubicCounter; [currentCubic, newCoeff]];
            else
                cubicCounter(matchingIdx,4) = cubicCounter(matchingIdx,4) + newCoeff;
            end
        end

        matchingIdx = repmat(vars,numQuartics,1);
        matchingIdx = matchingIdx == quarticCounter(:,1:4);
        matchingIdx = all(matchingIdx,2);
        matchingIdx = find(matchingIdx,1,'first');
        
        quarticCounter(matchingIdx,5) = quarticCounter(matchingIdx,5) + aj(vars(1))*aj(vars(2))*aj(vars(3)) * aj(vars(4)) / 16;
    end

end


matrixSize = size(x,1);
% form the singleton idx,
% these indices keep track which elements of basis x contain a given
% singleton
singletonGroupSizes = full(sum(x,1));
idxHelp = repmat(1:matrixSize,1,n);
singletonIdx = idxHelp( x(:) == 1);
singletonMonomials = (1:n)';
PRSMinput.singletonGroupSizes = singletonGroupSizes;
PRSMinput.singletonMonomials = singletonMonomials;
PRSMinput.singletonIdx = singletonIdx;
%%%% end singletonIdx creation

% compute all the cross products between monomials
% and transpose x for faster memory acces.
x = x';

numCrossPairs = 0.5*matrixSize*(matrixSize-1);
rowOffset = 1;

% store indices in vector allIdx
allIdx = zeros(numCrossPairs,1,'uint32');

% create simpleIdx to index the matrix used for constructing crossPairs
simpleIdx = allIdx;
simpleIdx(1) = 1;
for j = 2:matrixSize-1
    simpleIdx(rowOffset+1:rowOffset+j) = (1:j);
    rowOffset = rowOffset+j;
end

upperTriuIdx = simpleIdx + repelem(uint32(matrixSize:matrixSize:matrixSize*(matrixSize-1))',1:matrixSize-1);
% crossPairs is a huge matrix. Its orientation (i.e., with or without
% transpose) is carefully chosen as to optimize the memory acces speed
crossPairs = xor(repelem(x(:,2:end),1,1:matrixSize-1),x(:,simpleIdx));
monomialSizes = uint8(full(sum(crossPairs,1)));


groupSizes = zeros(numCrossPairs,1);
% groupSizes keeps track of the number of indices per monomial
if n > 4
    allMonomials = zeros(sum(monomialSizes),1);
else
    allMonomials = zeros(100,1);
end
varsCounter = 1;
idxCounter = 0;
groupSizeCounter = 1;
degCount = zeros(4,1);
% for k = 1


monSizeLogic = monomialSizes == 1;
currentPairs = crossPairs(:,monSizeLogic);
currentIdx = upperTriuIdx(monSizeLogic);
for p = 1:n
    matchingIdx = currentPairs(p,:);
    allMonomials(varsCounter) = p;
    varsCounter = varsCounter+1;

    idxCounter = idxCounter + 1;

    newIdx = currentIdx(matchingIdx);
    newGroupSize = size(newIdx,1);
    allIdx(groupSizeCounter:groupSizeCounter + newGroupSize -1) = newIdx;
    groupSizes(idxCounter) = newGroupSize;
    groupSizeCounter = groupSizeCounter + newGroupSize;
end
degCount(1) = n;
monomialSizes(monSizeLogic) = [];
crossPairs(:,monSizeLogic) = [];
upperTriuIdx(monSizeLogic) = [];

numCrossPairs = size(crossPairs,2);
groupUnitIdx = cumsum(groupSizes(1:idxCounter));
groupUnitIdx = groupUnitIdx(groupSizes(1:idxCounter) == 1);
unitCounter = size(groupUnitIdx,1)+1;
groupUnitIdx = [groupUnitIdx; zeros(numCrossPairs,1)];
groupIndicator = [repelem(uint32((1:idxCounter)'),groupSizes(1:idxCounter)); zeros(numCrossPairs,1)];
currentGroupNumber = idxCounter;

for k = 2:3
    monSizeLogic = monomialSizes == k;
    currentPairsBig = crossPairs(:,monSizeLogic);
    currentIdxBig = upperTriuIdx(monSizeLogic);
    for p = 1:(n-k+1)
        matchingVars = currentPairsBig(p,:);
        currentIdx  = currentIdxBig(matchingVars);
        if isempty(currentIdx)
            continue;
        end

        currentPairs = currentPairsBig(:,matchingVars);
        numCurrentPairs = size(currentPairs,2);

        % find the variables per clause
        allVars = find(currentPairs) - repelem(0:n:n*(numCurrentPairs-1),k)';
        allVars = reshape(allVars,k,numCurrentPairs)';
        allVars = allVars(:,2:k);

        % fast matlab 'unique' operator call
        [sortA,indSortA] = matlab.internal.math.sortrowsHelper(allVars);
        grpSortA = sortA(1:numCurrentPairs-1,:) ~= sortA(2:numCurrentPairs,:);
        grpSortA = [true; any(grpSortA,2)];
        uniqueVars = sortA(grpSortA,:);
        cIdx = indSortA(grpSortA);
        rIdx = cumsum(grpSortA);
        rIdx(indSortA) = rIdx;

        numUniqueRows = size(uniqueVars,1);
        uniqueVars = reshape([repmat(p,numUniqueRows,1), uniqueVars]',k*numUniqueRows,1);

        allMonomials(varsCounter:varsCounter+k*numUniqueRows-1) = uniqueVars;
        varsCounter = varsCounter+k*numUniqueRows;

        newGroupSizes = accumarray(rIdx,ones(numCurrentPairs,1));

        groupSizes(idxCounter+1:idxCounter+numUniqueRows) = newGroupSizes;
        idxCounter = idxCounter + numUniqueRows;
        unitGroups = cIdx(newGroupSizes == 1);
        if ~isempty(unitGroups)
            numNewUnits = size(unitGroups,1);
            groupUnitIdx(unitCounter:unitCounter+numNewUnits-1) = unitGroups+groupSizeCounter-1;
            unitCounter = unitCounter + numNewUnits;
        end

        groupIndicator(groupSizeCounter:groupSizeCounter+numCurrentPairs-1) = rIdx+currentGroupNumber;
        allIdx(groupSizeCounter:groupSizeCounter + numCurrentPairs-1) = currentIdx;
        currentGroupNumber = currentGroupNumber + numUniqueRows;
        groupSizeCounter = groupSizeCounter+ numCurrentPairs;

        currentPairsBig(:,matchingVars) = [];
        currentIdxBig(matchingVars) = [];
    end
    degCount(k) = idxCounter;
    
    monomialSizes(monSizeLogic) = [];
    crossPairs(:,monSizeLogic) = [];
    upperTriuIdx(monSizeLogic) = [];
end
% big matching loop end, k = 2:3



for k = 4
    % no need to first select only monomials with degree 4. All other
    % degrees have already been deleted from crossPairs and upperTriuIdx
    [monomialVector,~] = find(crossPairs);
    monomialVector = reshape(int16(monomialVector),k,size(crossPairs,2))';

    [~,~,ic] = unique(monomialVector,"rows");
    newIdx = upperTriuIdx;
    newGroupSizes = accumarray(ic,true(1,numel(ic)));
    
    unitGroups = newGroupSizes == 1;
    if any(unitGroups)
        unitGroups = find(unitGroups);
        unitGroupIdx = ismember(ic,unitGroups);
        newIdx(unitGroupIdx) = [];
        ic(unitGroupIdx) = [];
        newGroupSizes(unitGroups) = [];
    end
    groupSizes(idxCounter+1:idxCounter+numel(newGroupSizes)) = newGroupSizes;
    idxCounter = idxCounter + numel(newGroupSizes);

    groupIndicator(groupSizeCounter:groupSizeCounter+numel(newIdx)-1) = currentGroupNumber+uint32(ic);
    allIdx(groupSizeCounter:groupSizeCounter + numel(newIdx)-1) = newIdx;
    groupSizeCounter = groupSizeCounter+ numel(newIdx);


end
groupIndicator = groupIndicator(1:groupSizeCounter-1);
allIdx = allIdx(1:groupSizeCounter-1);
allMonomials = allMonomials(1:varsCounter-1);

if degCount(3)+1 <= degCount(4)
    idxVector = [1:degCount(1), repelem(degCount(1)+1:degCount(2),2),repelem(degCount(2)+1:degCount(3),3),repelem(degCount(3)+1:degCount(4),4)];
else
    idxVector = [1:degCount(1), repelem(degCount(1)+1:degCount(2),2),repelem(degCount(2)+1:degCount(3),3)];
end
allMonomials = sparse(allMonomials',idxVector ,true,n,idxCounter);
monomialsDegrees = repelem(uint8((1:3)'),[degCount(1), degCount(2)-degCount(1), degCount(3)-degCount(2)]);
groupSizes = groupSizes(1:idxCounter);
% only store the indices, not the coefficients since these are all 0
% anyways.
groupUnitIdx = groupUnitIdx(1:unitCounter-1);
% adept groupSizes, allIdx, groupIndicators

allIdx(groupUnitIdx) = [];
idxHelp = (1:idxCounter);

uniqueUnits = groupSizes == 1;
nonUnits = idxHelp(~uniqueUnits);
uniqueUnits = idxHelp(uniqueUnits);

groupIndicator(groupUnitIdx) = [];
% get the nonunit Idx monomials
% non unit, as in, they appear strictly more than 2 times in the symmetric
% matrix x*x'. Their indices are given by paddedNonUnitIdx
allNonUnitMon = allMonomials(:,nonUnits);



monomialsDegrees(uniqueUnits) = [];
groupSizes(uniqueUnits) = [];
% we match the coefficients of F_b to allMonomials
% match the terms corresponding to quadratics, i.e.,
% x_1*x_2, x_3*x_5. Determine here already to determine the size of coeffs
[row,col] = find(quadCounter);
quadIdx = [row,col]; quadCounter = quadCounter(quadCounter ~= 0) / 2;
coeffs = zeros(n+size(quadIdx,1)+size(cubicCounter,1)+size(quarticCounter,1),2);
% match the terms corresponding to singletons, i.e., x_1, x_2, x_3 etc
% divide by 2 for easier implementation
coeffs(1:n,:) = [(1:n)', linearCounter/2];

quadMonomials =  allNonUnitMon(:,monomialsDegrees == 2)';

numQuadMonomials = size(quadMonomials,1);
onesVector = true(2,1);

for j = 1:size(quadIdx,1)
    b = quadMonomials(:,quadIdx(j,:) ) * onesVector;
    matchingIdx = find(b == 2,1) + n;
    coeffs(n+j,:) = [matchingIdx, quadCounter(j)];
end
% now match the cubic terms to allMonomials
% divide by 2 for easier implementation
cubicCounter = sortrows(cubicCounter,[1 2 3]);
cubicCoeff = cubicCounter(:,4)/2;
cubicCounter = cubicCounter(:,1:3)';

numPrevMonomials = n + numQuadMonomials;
coeffIdxStart = n+size(quadIdx,1);

% match the degree 3 monomials to the right coefficient
% use that these are ordered
[g,~] = find(allNonUnitMon(:,monomialsDegrees == 3));
numCubicMonomials = numel(g)/3;
j = 1;
maxSize = size(cubicCounter,2);
if ~isempty(cubicCoeff)
    placeHolder = cubicCounter(:,j);
    for idx = 1:numCubicMonomials
        if all(g(3*idx-2:3*idx) == placeHolder)
            coeffs(coeffIdxStart+j,:) = [idx+numPrevMonomials, cubicCoeff(j)];
            j = j+1;
            if j > maxSize
                break;
            end
            placeHolder = cubicCounter(:,j);
        end
    end
end

% match the degree 4 monomials
quarticCounter = sortrows(quarticCounter,[1 2 3 4]);
quarticCoeff = quarticCounter(:,5)/2;
quarticCounter = quarticCounter(:,1:4)';

[g,~] = find(allNonUnitMon(:,monomialsDegrees == 4));
coeffIdxStart = coeffIdxStart + j-1;
j = 1;
maxSize = size(quarticCounter,2);
numPrevMonomials = numPrevMonomials+numCubicMonomials;

if ~isempty(quarticCoeff)
    placeHolder = quarticCounter(:,j);
    for idx = 1:size(g,1)/4
        if all(g(4*idx-3:4*idx) == placeHolder)
            coeffs(coeffIdxStart+j,:) = [idx+numPrevMonomials, quarticCoeff(j)];
            j = j+1;
            if j > maxSize
                break;
            end
            placeHolder = quarticCounter(:,j);
        end
    end
end



coeffs = sparse(coeffs(:,1),1,coeffs(:,2),size(allNonUnitMon,2),1);
% prepare output for PRSM function
PRSMinput.groupSizes = uint16(groupSizes);
PRSMinput.groupIndicators = uint32(groupIndicator);
PRSMinput.nonUnitIdx = allIdx;
PRSMinput.coeffs = coeffs;
PRSMinput.x = x';
PRSMinput.numClauses = m;
PRSMinput.numVariables = n;
PRSMinput.constant = constant;
end
