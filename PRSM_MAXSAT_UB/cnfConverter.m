function Q = cnfConverter(directoryName, fileName)
% Code accompanying the paper:
% On solving the MAX-SAT using sum of squares
% Lennart Sinjorgo, Renata Sotirov
% Feb 2023
% Tilburg, Netherlands
%
% Get Q matrix,
% Q has m rows, for every clause
% Q has 2n columns
% for i \in [n], Q(j,i) = 1 iff clause j contains variable i
%           nonnegated
% for i \in {n+1, n+2, ... , 2n}, Q(j,i) = 1 iff clause j contains variable
%           i negated.

currentDirectory = pwd;

% Change to directory containing CNF SAT instances
cd(directoryName);

% Enter fileName here
[fid,~] = fopen(fileName,'rt');
if fid == -1 
   cd(currentDirectory)
end
currentClause = 1;
locCounter = 1;
while true
  
  thisline = fgetl(fid);
  if ~ischar(thisline); break; end
  if thisline(1) == '%'; break; end
    if thisline(1) == 'p'
        % p cnf #variables #clauses
        % find the integers in the line
        num = regexp(thisline, '-?\d+', 'match');
        n = str2double(num{1});
        m = str2double(num{2});
        locations = zeros(round(n*m)/2,2);
    end
    if thisline(1) ~= 'c' && thisline(1) ~= 'p'
       % find the integers in the line, note that some integers
       % are negative. Last integer is always a 0, to indicate the end of
       % the line
       num = regexp(thisline, '-?\d+', 'match');
       if str2double(num{numel(num)}) ~= 0
          error('0 at the end of the line not found!') 
       end
       for i = 1:numel(num)-1
          int = str2double(num{i});
          if int > 0
             locations(locCounter,:) = [currentClause,abs(int)];
             locCounter = locCounter+1;
          end
          if int < 0
              locations(locCounter,:) = [currentClause,abs(int)+n];
              locCounter = locCounter+1;
          end
       end
       currentClause = currentClause+1;
    end
  end
  fclose(fid);
  cd(currentDirectory)
  
  locations = locations(1:locCounter-1,:);
  Q = sparse(locations(:,1),locations(:,2),true,m,2*n);
end
  