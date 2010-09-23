% diagonal scaling test
F = [ 1 0 ; 0 2 ];
G = [ 7 0 ; 0 8 ];
D = [ 5 0 ; 0 6 ];
C = [ 3 0 ; 0 4 ];

iF = inv(F);
iC = inv(C);

iAa = zeros(4,4);
iAa(1:2,1:2) = iF;
iAa(3:4,3:4) = iC;

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
C = [ 9  2 ;  6 5 ];

A = zeros(4,4);
A(1:2,1:2) = F;
A(1:2,3:4) = G;
A(3:4,1:2) = D;
A(3:4,3:4) = C;

iF = inv(F);
iC = inv(C);

iAa = zeros(4,4);
iAa(1:2,1:2) = iF;
iAa(3:4,3:4) = iC;

v = [0 1 1 3]';
o = (iAa*v)';
fprintf('\n   ef[0] = %.16e; ef[1] = %.16e;\n',o(1),o(2));
fprintf('   eg[0] = %.16e; eg[1] = %.16e;\n',o(3),o(4));

v = [-2 4 7 9]';
o = (iAa*v)';
fprintf('\n   ef[0] = %.16e; ef[1] = %.16e;\n',o(1),o(2));
fprintf('   eg[0] = %.16e; eg[1] = %.16e;\n',o(3),o(4));

v = [1 0 0 -5]';
o = (iAa*v)';
fprintf('\n   ef[0] = %.16e; ef[1] = %.16e;\n',o(1),o(2));
fprintf('   eg[0] = %.16e; eg[1] = %.16e;\n',o(3),o(4));

v = [4 -4 6 12]';
o = (iAa*v)';
fprintf('\n   ef[0] = %.16e; ef[1] = %.16e;\n',o(1),o(2));
fprintf('   eg[0] = %.16e; eg[1] = %.16e;\n',o(3),o(4));
