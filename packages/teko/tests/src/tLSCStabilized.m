% diagonal scaling test
F = [ 1 0 ; 0 2 ];
G = [ 7 0 ; 0 8 ];
D = [ 5 0 ; 0 6 ];
C = [ 3 0 ; 0 4 ];

gamma = 2.0/3.0;

N = D*inv(F)*G-C;
alpha = eig(D*inv(F)*G*inv(N));
alpha = 1.0/max(alpha);

iS = inv(D*G-gamma*C)*(D*F*G)*inv(D*G-gamma*C);
aiD = alpha*inv(N);

iU = eye(4);
iU(1:2,1:2) = inv(F);
iU(1:2,3:4) = inv(F)*G*(iS+aiD);
iU(3:4,3:4) = -iS-aiD;
iAa = iU;

v = [0 1 1 3]';
(iAa*v)'

v = [-2 4 7 9]';
(iAa*v)'

v = [1 0 0 -5]';
(iAa*v)'

v = [4 -4 6 12]';
(iAa*v)'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = diag([1,2,3,4,5]);
G = diag([7,8,9,10,11]);
D = diag([5,6,7,8,9]);
C = diag([3,4,5,6,7]);
Qu = diag([3,4,5,6,7]);
iQu = inv(Qu);

P = iQu*F/3.0;
gamma = max(P(:));
gamma
diag(inv(D*iQu*G-gamma*C),0)

BFBt = D*inv(F)*G;
M = BFBt-C;
iM = inv(M);
P = BFBt*iM;
alpha = 1.0/max(P(:))
diag(alpha*iM)

iS = inv(D*iQu*G-gamma*C)*D*iQu*F*iQu*G*inv(D*iQu*G-gamma*C) + alpha*iM;
iU = eye(10,10);
iU(1:5,1:5)   = inv(F);
iU(1:5,6:10)  = inv(F)*G*iS;
iU(6:10,6:10) = -iS;

iU
