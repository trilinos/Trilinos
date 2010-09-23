% build linear system to invert
A_00 = [1,2;2,1];
A_01 = [0,-1;3,4];
A_10 = [1,6;-3,2];
A_11 = [2,1;1,2];

% build vector
b_0 = [1,2];
b_1 = [3,4];
b = [b_0';b_1'];

% build pieces of preconditioner
alpha = 0.9;
P = A_00 + alpha*A_01;
H = diag(diag(A_11,0));

% build the preconditioner
M = [ P,0*A_01 ; A_10 , H];
iM = inv(M);

% compute action of preconditioner
x = iM*b
