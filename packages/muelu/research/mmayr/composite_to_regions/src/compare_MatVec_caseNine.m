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
%   |        |        |
%   |   r0   |   r1   |
%   |        |        |
%   +--------+--------+
%   |        |        |
%   |   r2   |   r3   |
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

% GIDs of each region
indReg0 = [1 2 3 6 7 8 11 12 13];
indReg1 = [3 4 5 8 9 10 13 14 15];
indReg2 = [11 12 13 16 17 18 21 22 23];
indReg3 = [13 14 15 18 19 20 23 24 25];

% LIDs of interface nodes ordered according to regions
lid0 = [3 6 7 8 9];
lid1 = [1 4 7 8 9];
lid2 = [1 2 3 6 9];
lid3 = [1 2 3 4 7];

% LIDs of non-Dirichlet nodes
nonDBC0 = [5 6 8 9];
nonDBC1 = [4 5 7 8];
nonDBC2 = [2 3 5 6];
nonDBC3 = [1 2 4 5];

% mapping of nodes to regions
nodesToRegion(1,:) = [1 -1 -1 -1];
nodesToRegion(2,:) = [1 -1 -1 -1];
nodesToRegion(3,:) = [1 2 -1 -1];
nodesToRegion(4,:) = [2 -1 -1 -1];
nodesToRegion(5,:) = [2 -1 -1 -1];
nodesToRegion(6,:) = [1 -1 -1 -1];
nodesToRegion(7,:) = [1 -1 -1 -1];
nodesToRegion(8,:) = [1 2 -1 -1];
nodesToRegion(9,:) = [2 -1 -1 -1];
nodesToRegion(10,:) = [2 -1 -1 -1];
nodesToRegion(11,:) = [1 3 -1 -1];
nodesToRegion(12,:) = [1 3 -1 -1];
nodesToRegion(13,:) = [1 2 3 4];
nodesToRegion(14,:) = [2 4 -1 -1];
nodesToRegion(15,:) = [2 4 -1 -1];
nodesToRegion(16,:) = [3 -1 -1 -1];
nodesToRegion(17,:) = [3 -1 -1 -1];
nodesToRegion(18,:) = [3 4 -1 -1];
nodesToRegion(19,:) = [4 -1 -1 -1];
nodesToRegion(20,:) = [4 -1 -1 -1];
nodesToRegion(21,:) = [3 -1 -1 -1];
nodesToRegion(22,:) = [3 -1 -1 -1];
nodesToRegion(23,:) = [3 4 -1 -1];
nodesToRegion(24,:) = [4 -1 -1 -1];
nodesToRegion(25,:) = [4 -1 -1 -1];

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
regA3 = scaledA(indReg3, indReg3);

%% Fix diagonal values

if (strcmp(nullspaceScheme, 'exact preservation') == true)
  % Modify diagonal entries to preserve nullspace. Since we want to preserve the
  % nullspace in the region layout, this is done region by region.
  regX0 = ones(9,1); % nullspace vector
  tmp = regA0 * regX0; % compute action of matrix on nullspace vector
  for i = nonDBC0 % only process non-DBC nodes
    regA0(i,i) = regA0(i,i) - tmp(i); % correct diagonal entry
  end

  % compute violation of nullspace preservation and correct
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

  regX3 = ones(9,1);
  tmp = regA3 * regX3;
  for i = nonDBC3
    regA3(i,i) = regA3(i,i) - tmp(i);
  end

  % Process interface nodes subject to Dirichlet boundary conditions
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
  for i = lid3
    if regA3(i,i) == 1
      regA3(i,i) = 0.5;
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
  nspViolationReg3 = compViolation(indReg3);
  
  nspViolationReg0(lid0) = nspViolationReg0(lid0) / 2;
  nspViolationReg1(lid1) = nspViolationReg1(lid1) / 2;
  nspViolationReg2(lid2) = nspViolationReg2(lid2) / 2;
  nspViolationReg3(lid3) = nspViolationReg3(lid3) / 2;
  
  regX0 = ones(9,1); % nullspace vector
  tmp = regA0 * regX0; % compute action of matrix on nullspace vector
  for i = 1:9 % process all diagonal entries
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
  
  regX3 = ones(9,1);
  tmp = regA3 * regX3;
  for i = 1:9
    regA3(i,i) = regA3(i,i) - tmp(i) + nspViolationReg3(i);
  end
  
else
  error('Unknown nullspace scheme "%s".', nullspaceScheme);
end

%% Build global, but regional matrices and vectors
% assemble regA from regional submatrices
regZS = zeros(9,9);
regA = [regA0 regZS regZS regZS;
        regZS regA1 regZS regZS;
        regZS regZS regA2 regZS;
        regZS regZS regZS regA3];

% extract regional vectors from composite vector and assemle vector in region
% layout
regX0 = x(indReg0);
regX1 = x(indReg1);
regX2 = x(indReg2);
regX3 = x(indReg3);
regX = [regX0; regX1; regX2; regX3];

%% Perform matrix-vector multiplication in region layout
regZ = regA * regX;

%% Transform result back to composite regime and compare to composite result

regZ0 = regZ(1:9);
regZ1 = regZ(10:18);
regZ2 = regZ(19:27);
regZ3 = regZ(28:36);

scaledCompZ = zeros(size(z));

scaledCompZ(1) = regZ0(1);
scaledCompZ(2) = regZ0(2);
scaledCompZ(3) = regZ0(3) + regZ1(1);
scaledCompZ(4) = regZ1(2);
scaledCompZ(5) = regZ1(3);
scaledCompZ(6) = regZ0(4);
scaledCompZ(7) = regZ0(5);
scaledCompZ(8) = regZ0(6) + regZ1(4);
scaledCompZ(9) = regZ1(5);
scaledCompZ(10) = regZ1(6);
scaledCompZ(11) = regZ0(7) + regZ2(1);
scaledCompZ(12) = regZ0(8) + regZ2(2);
scaledCompZ(13) = regZ0(9) + regZ1(7) + regZ2(3) + regZ3(1);
scaledCompZ(14) = regZ1(8) + regZ3(2);
scaledCompZ(15) = regZ1(9) + regZ3(3);
scaledCompZ(16) = regZ2(4);
scaledCompZ(17) = regZ2(5);
scaledCompZ(18) = regZ2(6) + regZ3(4);
scaledCompZ(19) = regZ3(5);
scaledCompZ(20) = regZ3(6);
scaledCompZ(21) = regZ2(7);
scaledCompZ(22) = regZ2(8);
scaledCompZ(23) = regZ2(9) + regZ3(7);
scaledCompZ(24) = regZ3(8);
scaledCompZ(25) = regZ3(9);

str = sprintf('GID\tz\tscaledCompZ\tdiff');
disp(str);
disp([[1:25]' z scaledCompZ z-scaledCompZ]);

str = sprintf('norm of difference: %f', norm(z-scaledCompZ));
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