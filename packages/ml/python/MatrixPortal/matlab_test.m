n = 10;

e = ones(n,1);
A = spdiags([e -2*e e], -1:1, n, n);

ExactSolution = rand(n, 1);
RHS = A * ExactSolution;
LHS = zeros(n, 1);
NullSpace = rand(n, 1); % bad choice, only for testing

matlab_hdf5('matlab_test.h5', A, RHS, LHS, ExactSolution, NullSpace);
