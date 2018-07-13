function [X, nnzlu] = LUdriver(A, k, B)
%function X = LUdriver(A, k, B)
%
% A driver for the LU factorization using hypergraph partitioner.
%
% Input:
%   A, matrix to be factored.
%   k, (integer) number of parts.
%   B, right hand side
% Output:
%   X, left hand side
%   nnzlu, not including the non zeros in the band of A, including the non
%   zeros in the factor of the Schur complement.
%
% Erik Boman, Sandia National Labs, 2010.
% Siva Rajamanickam, Sandia National Labs, 2010.

if (nargin<2)
  k = max(2, ceil(size(A,1)/1000)); % Default value for #blocks
end

[m,n] = size(A);

if (m ~= n)
  error('Input A must be square');
  return;
end

if (k<1)
  error('Input k must be > 1');
  return;
end

if (nargin < 3)
  B = sprand(n, 1, 0.4);
end

% Compute hypergraph partition and the LU factors of the diagonal blocks.
LUdata = LUfactor(A, k);
% Find Bhat, the rhs for Schur complement.
Bhat = LUfsolve(LUdata, B);
nnzlu = LUdata.nnzlu ;

khalf = floor(k/2) ;

if (khalf == 1)
  % Solve Schur complement. Sequential bottle neck.
  LUdata.Xhat = LUdata.U{k+1} \ (LUdata.L{k+1} \ Bhat(LUdata.p{k+1}));

  % Solve the diagonal blocks
  X = LUbsolve(LUdata, B);
else
  % Ignore the L and U of the Schur complement, refactor S using hypergraph
  % partitioning
    %%SLUdata = LUfactor(LUdata.S{k+1}, khalf);
    %%B1hat = LUfsolve(SLUdata, Bhat);
    %%SLUdata.Xhat = SLUdata.U{khalf+1} \ (SLUdata.L{khalf+1} \ ...
										%%B1hat(SLUdata.p{khalf+1}));
    %%LUdata.Xhat = LUbsolve(SLUdata, Bhat);

  LUdata.q{k+1} = LUdata.c{k+1};
  [LUdata.Xhat, nnzlu1] = LUdriver(LUdata.S{k+1}, khalf, Bhat) ;
  % Adjust for the non zeros the factors of S{k+1}
  nnzlu = nnzlu - nnz(LUdata.L{k+1}) - nnz(LUdata.U{k+1}) + nnzlu1 ;
  % Solve the diagonal blocks.
  X = LUbsolve(LUdata, B);
end
