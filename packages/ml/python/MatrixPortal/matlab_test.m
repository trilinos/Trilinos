n = 10;

e = ones(n,1);
A = spdiags([e -2*e e], -1:1, n, n)

ExactSolution = rand(n, 1);
RHS = A * ExactSolution;
LHS = zeros(n, 1);
NullSpace = ones(n, 1);

matlab_hdf5('matlab.h5', A);
