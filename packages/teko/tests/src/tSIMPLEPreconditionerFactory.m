% diagonal scaling test
F = [ 1 0 ; 0 2 ];
G = [ 7 0 ; 0 8 ];
D = [ 5 0 ; 0 6 ];
C = [ 3 0 ; 0 4 ];

alpha = 0.9;
H = diag(1./diag(F),0);
S = -C + D*H*G

L = [ F, 0*G ; D, -S ];
U = [ eye(size(F)), H*G/alpha ; 0*D, eye(size(C))/alpha ];
A_S = L*U;
iA = inv(A_S);

v = [0 1 1 3]';
(iA*v)'

v = [-2 4 7 9]';
(iA*v)'

v = [1 0 0 -5]';
(iA*v)'

v = [4 -4 6 12]';
(iA*v)'

disp('*************************************');
clear;

% fully defined system test
F = [ 1  2 ;  2 1 ];
G = [ 1 -1 ; -3 1 ];
D = [ 1 -3 ; -1 1 ];
C = [ 1  2 ;  2 1 ];

alpha = 0.9;
H = diag(1./diag(F),0);
S = -C + D*H*G;

inv(S)

L = [ F, 0*G ; D, -S ];
U = [ eye(size(F)), H*G/alpha ; 0*D, eye(size(C))/alpha ];
A_S = L*U;
iA = inv(A_S);

v = [0 1 1 3]';
(iA*v)'

v = [-2 4 7 9]';
(iA*v)'

v = [1 0 0 -5]';
(iA*v)'

v = [4 -4 6 12]';
(iA*v)'
