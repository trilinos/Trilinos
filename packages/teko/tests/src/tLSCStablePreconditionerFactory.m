% diagonal scaling test
F = [ 1 0 ; 0 2 ];
G = [ 7 0 ; 0 8 ];
D = [ 5 0 ; 0 6 ];

iS = inv(D*G)*(D*F*G)*inv(D*G);

%iL = eye(4);
%iD = zeros(4);
%iU = eye(4);

%iD(1:2,1:2) = inv(F);
%iD(3:4,3:4) = -iS;

%iU(1:2,3:4) = -inv(F)*G;
%iL(3:4,1:2) = -D*inv(F);
%
%iAa = iU*iD*iL;

iU = eye(4);
iU(1:2,1:2) = inv(F);
iU(1:2,3:4) = inv(F)*G*iS;
iU(3:4,3:4) = -iS;
iAa = iU;

v = [0 1 1 3]';
(iAa*v)'

v = [-2 4 7 9]';
(iAa*v)'

v = [1 0 0 -5]';
(iAa*v)'

v = [4 -4 6 12]';
(iAa*v)'

disp('*************************************');
clear;

% fully defined system test
F = [ 1  2 ;  2 1 ];
G = [ 1 -1 ; -3 1 ];
D = [ 1 -3 ; -1 1 ];
iM = [ 3 0 ; 0 2 ];

A = zeros(4,4);
A(1:2,1:2) = F;
A(1:2,3:4) = G;
A(3:4,1:2) = D;

iBQBt = inv(D*iM*G);
iS = iBQBt*(D*iM*F*iM*G)*iBQBt;

% iL = eye(4);
% iD = zeros(4);
% iU = eye(4);
% 
% iD(1:2,1:2) = inv(F);
% iD(3:4,3:4) = -iS;
% 
% iU(1:2,3:4) = -inv(F)*G;
% iL(3:4,1:2) = -D*inv(F);
% 
% iAa = iU*iD*iL;

iU = eye(4);
iU(1:2,1:2) = inv(F);
iU(1:2,3:4) = inv(F)*G*iS;
iU(3:4,3:4) = -iS;
iAa = iU;

v = [0 1 1 3]';
(iAa*v)'

v = [-2 4 7 9]';
(iAa*v)'

v = [1 0 0 -5]';
(iAa*v)'

v = [4 -4 6 12]';
(iAa*v)'
