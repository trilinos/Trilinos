% Compute matrix-vector product for caseNine.
%
% This mimics a region-wise matrix-vector product for a two-dimensional problem.
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

%% Define region-wise quantities as extraction from global matrix

% region 0
regA0 = eye(9);
regA0(5,:) = [-1 -1 -1 -1  8 -1 -1 -1 -1];
regA0(6,:) = [ 0 -1 -1  0 -1  8  0 -1 -1];
regA0(8,:) = [ 0  0  0 -1 -1 -1 -1  8 -1];
regA0(9,:) = [ 0  0  0  0 -1 -1  0 -1  8];

% region 1
regA1 = eye(9);
regA1(4,:) = [-1 -1  0  8 -1  0 -1 -1  0];
regA1(5,:) = [-1 -1 -1 -1  8 -1 -1 -1 -1];
regA1(7,:) = [ 0  0  0 -1 -1  0  8 -1  0];
regA1(8,:) = [ 0  0  0 -1 -1 -1 -1  8 -1];

% region 2
regA2 = eye(9);
regA2(2,:) = [-1  8 -1 -1 -1 -1  0  0  0];
regA2(3,:) = [ 0 -1  8  0 -1 -1  0  0  0];
regA2(5,:) = [-1 -1 -1 -1  8 -1 -1 -1 -1];
regA2(6,:) = [ 0 -1 -1  0 -1  8  0 -1 -1];

% region 3
regA3 = eye(9);
regA3(1,:) = [ 8 -1  0 -1 -1  0  0  0  0];
regA3(2,:) = [-1  8 -1 -1 -1 -1  0  0  0];
regA3(4,:) = [-1 -1  0  8 -1  0 -1 -1  0];
regA3(5,:) = [-1 -1 -1 -1  8 -1 -1 -1 -1];

%% Modify entries in interface matrices

% % local IDs of interface nodes per region
% int0 = [3 6 7 8 9];
% int1 = [1 4 7 8 9];
% int2 = [1 2 3 6 9];
% int3 = [1 2 3 4 7];

% for i = int0
%   regA0(i,int0) = regA0(i,int0) / 2;
% end
% for i = int1
%   regA1(i,int1) = regA1(i,int1) / 2;
% end
% for i = int2
%   regA2(i,int2) = regA2(i,int2) / 2;
% end
% for i = int3
%   regA3(i,int3) = regA3(i,int3) / 2;
% end
% 
% 
% % fix center node that has been duplicated three times
% regA0(9,9) = regA0(9,9) * 2 / 4;
% regA1(7,7) = regA1(7,7) * 2 / 4;
% regA2(3,3) = regA2(3,3) * 2 / 4;
% regA3(1,1) = regA3(1,1) * 2 / 4;

% local IDs of interface nodes per region and interface ID
int00 = [3 6 9];
int01 = [7 8 9];
int10 = [1 4 7];
int12 = [7 8 9];
int21 = [1 2 3];
int23 = [3 6 9];
int32 = [1 2 3];
int33 = [1 4 7];

for i = int00
  regA0(i, int00) = regA0(i, int00) / 2;
end
for i = int01
  regA0(i, int01) = regA0(i, int01) / 2;
end
for i = int10
  regA1(i, int10) = regA1(i, int10) / 2;
end
for i = int12
  regA1(i, int12) = regA1(i, int12) / 2;
end
for i = int21
  regA2(i, int21) = regA2(i, int21) / 2;
end
for i = int23
  regA2(i, int23) = regA2(i, int23) / 2;
end
for i = int32
  regA3(i, int32) = regA3(i, int32) / 2;
end
for i = int33
  regA3(i, int33) = regA3(i, int33) / 2;
end


%% Build global, but regional matrices

regZS = zeros(9,9);

regA = [regA0 regZS regZS regZS;
        regZS regA1 regZS regZS;
        regZS regZS regA2 regZS;
        regZS regZS regZS regA3];
regX = ones(size(regA,1),1);

%% composite quantities
A = full(twoDimensionalLaplace(25));

x = ones(size(A,1),1);


%% Perform matrix-vector multiplication
z = A*x;

regZ = regA * regX;

%% Transform result back to composite regime

regZ0 = regZ(1:9);
regZ1 = regZ(10:18);
regZ2 = regZ(19:27);
regZ3 = regZ(28:36);

compZ = zeros(size(z));

compZ(1) = regZ0(1);
compZ(2) = regZ0(2);
compZ(3) = regZ0(3) + regZ1(1);
compZ(4) = regZ1(2);
compZ(5) = regZ1(3);
compZ(6) = regZ0(4);
compZ(7) = regZ0(5);
compZ(8) = regZ0(6) + regZ1(4);
compZ(9) = regZ1(5);
compZ(10) = regZ1(6);
compZ(11) = regZ0(7) + regZ2(1);
compZ(12) = regZ0(8) + regZ2(2);
compZ(13) = regZ0(9) + regZ1(7) + regZ2(3) + regZ3(1);
compZ(14) = regZ1(8) + regZ3(2);
compZ(15) = regZ1(9) + regZ3(3);
compZ(16) = regZ2(4);
compZ(17) = regZ2(5);
compZ(18) = regZ2(6) + regZ3(4);
compZ(19) = regZ3(5);
compZ(20) = regZ3(6);
compZ(21) = regZ2(7);
compZ(22) = regZ2(8);
compZ(23) = regZ2(9) + regZ3(7);
compZ(24) = regZ3(8);
compZ(25) = regZ3(9);

str = sprintf('GID\tz\tcomZ');
disp(str);
disp([[1:25]' z compZ]);

str = sprintf('norm of difference: %f', norm(z-compZ));
disp(str);
