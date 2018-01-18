% Exercise edge-based splitting with nullspace constraint
%
% We consider an edge-based approach. We assume that each matrix entry a_ij 
% (i != j) represents an 'edge' that connects the nodes i and j. It is known, to 
% which regions i and j belong to, respsectiely. 
% - If i and j belong to one and the same region, this is an interior edge.
% - If i and j belong to the same two regions, this is an interface edge.
%
% Diagonal entries a_ii are not considered to be edges, but require a differnet
% treatment.
%
% The splitting algorithm consists of two steps:
% 1. Process off-diagonals: Interior edges don't need to be changed. Interface 
%    edges are shared by 2 regions and, thus, need to be splitted (divided by 2).
% 2. Process diagonal entries: Diagonal entries need to be chosen to preserve
%    nullspace properties for region matrices. We take the kth region matrix
%    regAk and multiply it with its portion of the global nullspace vector, i.e.
%    tmp = regAk * regNSPk. The result tmp now contains the correction for the
%    diagonal entries, i.e. diag(regAk) = diag(regAk) - tmp.
%
% Region numbering:
%
%   +--------+--------+
%   |                 |
%   |       r0        |
%   |                 |
%   +--------+--------+
%   |        |        |
%   |   r1   |   r2   |
%   |        |        |
%   +--------+--------+
%

clear;

%% User's choices

% nullspaceScheme = 'exact preservation';
nullspaceScheme = 'preserve violation';

%% composite quantities
A = full(twoDimensionalLaplace(25));

% x = ones(size(A,1),1);
% x = [1:25]';
x = rand(25,1);

z = A * x;


%% Define region-wise quantities as extraction from global matrix

% work on a copy of the composite matrix for scaling
scaledA = A;

indReg0 = [1:15];
indReg1 = [11 12 13 16 17 18 21 22 23];
indReg2 = [13 14 15 18 19 20 23 24 25];

% LIDs of interface nodes ordered according to regions
lid0 = [11 12 13 14 15];
lid1 = [1 2 3 6 9];
lid2 = [1 2 3 4 7];

% LIDs of non-Dirichlet nodes
nonDBC0 = [7 8 9 12 13 14];
nonDBC1 = [2 3 5 6];
nonDBC2 = [1 2 4 5];

% mapping of nodes to regions
nodesToRegion(1,:)  = [1 -1 -1];
nodesToRegion(2,:)  = [1 -1 -1];
nodesToRegion(3,:)  = [1 -1 -1];
nodesToRegion(4,:)  = [1 -1 -1];
nodesToRegion(5,:)  = [1 -1 -1];
nodesToRegion(6,:)  = [1 -1 -1];
nodesToRegion(7,:)  = [1 -1 -1];
nodesToRegion(8,:)  = [1 -1 -1];
nodesToRegion(9,:)  = [1 -1 -1];
nodesToRegion(10,:) = [1 -1 -1];
nodesToRegion(11,:) = [1  2 -1];
nodesToRegion(12,:) = [1  2 -1];
nodesToRegion(13,:) = [1  2  3];
nodesToRegion(14,:) = [1  3 -1];
nodesToRegion(15,:) = [1  3 -1];
nodesToRegion(16,:) = [2 -1 -1];
nodesToRegion(17,:) = [2 -1 -1];
nodesToRegion(18,:) = [2  3 -1];
nodesToRegion(19,:) = [3 -1 -1];
nodesToRegion(20,:) = [3 -1 -1];
nodesToRegion(21,:) = [2 -1 -1];
nodesToRegion(22,:) = [2 -1 -1];
nodesToRegion(23,:) = [2  3 -1];
nodesToRegion(24,:) = [3 -1 -1];
nodesToRegion(25,:) = [3 -1 -1];

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
  regX0 = ones(15,1); % nullspace vector
  tmp = regA0 * regX0; % compute action of matrix on nullspace vector
  for i = nonDBC0 % only process non-DBC nodes
    regA0(i,i) = regA0(i,i) - tmp(i); % correct diagonal entry
  end

  regX1 = ones(9,1);
  tmp = regA1 * regX1;
  for i = nonDBC1
    regA1(i,i) = regA1(i,i) - tmp(i);
  end

  regX2 = ones(9,1);
  tmp = regA2 * regX2;
  for i = nonDBC2
    regA2(i,i) = regA2(i,i) - tmp(i);
  end

  for i = lid0
    if regA0(i,i) == 1
      regA0(i,i) = 0.5;
    end
  end
  for i = lid1
    if regA1(i,i) == 1
      regA1(i,i) = 0.5;
    end
  end
  for i = lid2
    if regA2(i,i) == 1
      regA2(i,i) = 0.5;
    end
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
  
  regX0 = ones(15,1); % nullspace vector
  tmp = regA0 * regX0; % compute action of matrix on nullspace vector
  for i = 1:15 % process all diagonal entries
    regA0(i,i) = regA0(i,i) - tmp(i) + nspViolationReg0(i); % correct diagonal entry
  end
    
  regX1 = ones(9,1);
  tmp = regA1 * regX1;
  for i = 1:9
    regA1(i,i) = regA1(i,i) - tmp(i) + nspViolationReg1(i);
  end
  
  regX2 = ones(9,1);
  tmp = regA2 * regX2;
  for i = 1:9
    regA2(i,i) = regA2(i,i) - tmp(i) + nspViolationReg2(i);
  end
   
else
  error('Unknown nullspace scheme "%s".', nullspaceScheme);
end

%% Build global, but regional matrices and vectors
% assemble regA from regional submatrices
regZ00 = zeros(15,15);
regZ01 = zeros(15,9);
regZ10 = regZ01';
regZ11 = zeros(9,9);

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

regZ0 = regZ(1:15);
regZ1 = regZ(16:24);
regZ2 = regZ(25:33);

compZ = zeros(size(z));

compZ(1:10) = regZ0(1:10);
compZ(11) = regZ0(11) + regZ1(1);
compZ(12) = regZ0(12) + regZ1(2);
compZ(13) = regZ0(13) + regZ1(3) + regZ2(1);
compZ(14) = regZ0(14) + regZ2(2);
compZ(15) = regZ0(15) + regZ2(3);
compZ([16 17 21 22]) = regZ1([4 5 7 8]);
compZ([19 20 24 25]) = regZ2([5 6 8 9]);
compZ([18 23]) = regZ1([6 9]) + regZ2([4 7]);

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