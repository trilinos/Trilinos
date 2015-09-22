n = 10;
N = n*n;

% Create 2D Laplacian matrix
I = speye(n,n);
E = sparse(2:n, 1:n-1, 1, n, n);
D = E+E'- 2*I;
A = kron(D,I) + kron(I,D);

% Create rhs and initial guess
b = [N:-1:1]';
x = zeros(N,1);

% Do the setup stage
[h, oc] = muelu('setup', A);
