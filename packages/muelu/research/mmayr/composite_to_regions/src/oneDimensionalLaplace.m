% oneDimensionalLaplace.m
%
% Build FE stiffness matrix for 1D Laplace problem with Dirichlet BC on the left
% and Neumann BC on the right.
%
function [ A ] = oneDimensionalLaplace(nNodes)

A = zeros(nNodes,nNodes);
A(1,1) = 1.0;
for i = 2:nNodes-1
  A(i,i-1) = -1.0;
  A(i,i) = 2.0;
  A(i,i+1) = -1.0;
end
A(end,end - 1) = -1.0;
A(end,end) = 1.0;

A(1,2) = -1.0e-14;

A = sparse(A);

end