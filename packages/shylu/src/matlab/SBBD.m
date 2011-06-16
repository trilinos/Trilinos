function [partr, partc, permr, permc, domain] = SBBD(A, k)
%function [partr, partc, permr, permc, domain] = SBBD(A, k)
%
% Partition/permute a sparse matrix to
% singly bordered block diagonal (SBBD) form.
%
% Note this function needs PaToH to get good results!
%
% Input:
%   A, a sparse matrix to partition/permute
%   k, the desired number of diagonal blocks
%
% Output:
%   partr: part numbers for rows (1:k)
%   partc: part numbers for columns (1:k+1)
%   permr: permutation vector for rows
%   permc: permutation vector for columns
%
% Author: Erik Boman, Sandia National Labs, 2010.

if (k < 1)
  error('Input argument k must be at least 1.');
  return;
end

[m,n] = size(A);

if exist('PaToH')
  % Use Patoh to partition rows
  partr = PaToH(A, k) +1; % Make 1-based
else
  % Simple block partitioning
  partr = ceil(k*(1:m)/m);
end

% Support older version of Patoh that returns row vectors.
[nr, nc] = size(partr);
if (nc ~= 1)
  partr = partr';
end

% Compute column blocks
% Modify A so we can find internal and cut nets (columns)
A = spones(A);
% TODO: Row scaling as matrix multiply?
for i= 1:m
  A(i,:) = partr(i) * A(i,:);
end

% Find min and max part for each column
maxpart = max(A);
%minpart = min(A);% This gives 0; we want the min of the nonzeros
minpart= zeros(size(maxpart));
[ii,jj,val] = find(A);
for j= 1:n
  minpart(j) = min(val(find(jj==j)));
end
partc = zeros(n,1);
domain = zeros(k, k);
for j= 1:n
  x = zeros(k,1);
  if (minpart(j)==maxpart(j))
    partc(j) = minpart(j);
  else
    partc(j) = k+1; % Cut net (column)
    dno = unique(val(find(jj == j))) ;
    x(dno) = 1;
    domain = domain + (x * x') ;
  end
end

% If requested, return permutations
if (nargout > 2)
  [dummy, permr] = sort(partr);
  [dummy, permc] = sort(partc);
end
