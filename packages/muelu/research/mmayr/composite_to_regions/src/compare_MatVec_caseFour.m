% Compute coarse level operator for caseFour.
%
% Exercise the edge-based splitting approach for one-dimensional problems

clear;

%% User's choices

nullspaceScheme = 'exact preservation';
% nullspaceScheme = 'preserve violation';

%% composite quantities
A = full(oneDimensionalLaplace(25));

% x = ones(size(A,1),1);
% x = [1:25]';
x = rand(25,1);

z = A * x;

%% Define region-wise quantities as extraction from global matrix

% work on a copy of the composite matrix for scaling
scaledA = A;

% GIDs of each region
indReg0 = [1:7];
indReg1 = [7:16];
indReg2 = [16:25];

% LIDs of interface nodes ordered according to regions
lid0 = [7];
lid1 = [1 10];
lid2 = [1];

% LIDs of non-Dirichlet nodes
nonDBC0 = [2:7];
nonDBC1 = [1:10];
nonDBC2 = [1:10];

% mapping of nodes to regions
nodesToRegion = [1 -1;
                 1 -1;
                 1 -1;
                 1 -1;
                 1 -1;
                 1 -1;
                 1  2;
                 2 -1;
                 2 -1;
                 2 -1;
                 2 -1;
                 2 -1;
                 2 -1;
                 2 -1;
                 2 -1;
                 2  3;
                 3 -1;
                 3 -1;
                 3 -1;
                 3 -1;
                 3 -1;
                 3 -1;
                 3 -1;
                 3 -1;
                 3 -1];
               
% We first process the off-diagonals. For convenience of GID access, this is
% done in the composite matrix.
for r = 1:size(scaledA,1) % loop over all rows
  row = scaledA(r,:); % grep a single row
  [~, nnzInds, ~] = find(row); % identify nonzeros in this row
  for c = nnzInds % loop over all nonzeros in this row
    if r ~= c % skip the diagonal element
      commonRegs = findCommonRegions(nodesToRegion(r,:), nodesToRegion(c,:));
      numCommonRegs = length(commonRegs);
      
      if numCommonRegs == 2
        scaledA(r,c) = scaledA(r,c) / 2;
      elseif numCommonRegs == 0 || numCommonRegs == 1
        % do nothing
      else
        error('Search for common region assignments produced weird result.');
      end
    end
  end
end

% extract region matrices
regA0 = scaledA(indReg0, indReg0);
regA1 = scaledA(indReg1, indReg1);
regA2 = scaledA(indReg2, indReg2);

%% Fix diagonal values

if (strcmp(nullspaceScheme, 'exact preservation') == true)
  % Modify diagonal entries to preserve nullspace. Since we want to preserve the
  % nullspace in the region layout, this is done region by region.
  regX0 = ones(7,1); % nullspace vector
  tmp = regA0 * regX0; % compute action of matrix on nullspace vector
  for i = nonDBC0 % only process non-DBC nodes
    regA0(i,i) = regA0(i,i) - tmp(i); % correct diagonal entry
  end

  regX1 = ones(10,1);
  tmp = regA1 * regX1;
  for i = nonDBC1
    regA1(i,i) = regA1(i,i) - tmp(i);
  end

  regX2 = ones(10,1);
  tmp = regA2 * regX2;
  for i = nonDBC2
    regA2(i,i) = regA2(i,i) - tmp(i);
  end
  
elseif (strcmp(nullspaceScheme, 'preserve violation') == true)
  % compute violation of nullspace preservation in the composite layout
  compNsp = ones(size(x));
  compViolation = A * compNsp;
  
  % extract region-wise violation
  nspViolationReg0 = compViolation(indReg0);
  nspViolationReg1 = compViolation(indReg1);
  nspViolationReg2 = compViolation(indReg2);
    
  nspViolationReg0(lid0) = nspViolationReg0(lid0) / 2;
  nspViolationReg1(lid1) = nspViolationReg1(lid1) / 2;
  nspViolationReg2(lid2) = nspViolationReg2(lid2) / 2;
    
  regX0 = ones(7,1); % nullspace vector
  tmp = regA0 * regX0; % compute action of matrix on nullspace vector
  for i = 1:7 % process all diagonal entries
    regA0(i,i) = regA0(i,i) - tmp(i) + nspViolationReg0(i); % correct diagonal entry
  end
    
  regX1 = ones(10,1);
  tmp = regA1 * regX1;
  for i = 1:10
    regA1(i,i) = regA1(i,i) - tmp(i) + nspViolationReg1(i);
  end
  
  regX2 = ones(10,1);
  tmp = regA2 * regX2;
  for i = 1:10
    regA2(i,i) = regA2(i,i) - tmp(i) + nspViolationReg2(i);
  end
   
else
  error('Unknown nullspace scheme "%s".', nullspaceScheme);
end

%% Build global, but regional matrices and vectors
% assemble regA from regional submatrices
regZ00 = zeros(7,7);
regZ01 = zeros(7,10);
regZ10 = regZ01';
regZ11 = zeros(10,10);

regA = [regA0 regZ01 regZ01;
        regZ10 regA1 regZ11;
        regZ10 regZ11 regA2];
      
% extract regional vectors from composite vector and assemle vector in region
% layout
regX0 = x(indReg0);
regX1 = x(indReg1);
regX2 = x(indReg2);
regX = [regX0; regX1; regX2];

%% Perform matrix-vector multiplication in region layout
regZ = regA * regX;

%% Transform result back to composite regime and compare to composite result

regZ0 = regZ(1:7);
regZ1 = regZ(8:17);
regZ2 = regZ(18:27);

compZ = zeros(size(z));

compZ(1:6) = regZ0(1:end-1);
compZ(7) = regZ0(end) + regZ1(1);
compZ(8:15) = regZ1(2:end-1);
compZ(16) = regZ1(end) + regZ2(1);
compZ(17:25) = regZ2(2:end);

str = sprintf('GID\tz\tcompZ');
disp(str);
disp([[1:25]' z compZ]);

str = sprintf('norm of difference: %f', norm(z-compZ));
disp(str);

%% Utils: find shared regions of two nodes
function [ commonRegs ] = findCommonRegions(regsA, regsB)

commonRegs = [];

for a = regsA
  for b = regsB
    if (a == b)
      commonRegs = [commonRegs a];
    end
  end
end

commonRegs = commonRegs(find(commonRegs > 0));
commonRegs = unique(commonRegs);

end
