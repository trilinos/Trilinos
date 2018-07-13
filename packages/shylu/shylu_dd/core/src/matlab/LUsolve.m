function X = LUsolve(LUdata, B)
% function X = LUsolve(LUdata, B)
%
% 2-level LU based on SBBD form.
% Solve phase only. LUfactor must be called first.
%
% Erik Boman, Sandia National Labs, 2010.

[n,nrhs] = size(B);
k = LUdata.nblocks;
A = LUdata.A;

% Modify rhs for Schur complement.
Bhat = [];
for i= 1:k
  % b_i = b_i - \bar{A_i} * (A_i \ b_i)
  Bi = B(LUdata.Sr{i},:) - (A(LUdata.Sr{i},LUdata.q{i})/LUdata.U{i}) * ...
		(LUdata.L{i}\ B(LUdata.p{i},:));
  Bhat = [Bhat; Bi];
end

% Solve Schur complement.
Xhat{k+1} = LUdata.U{k+1} \ (LUdata.L{k+1} \ Bhat(LUdata.p{k+1}));

% Solve for diagonal blocks.
for i= 1:k
  % b_i = b_i - A_{i,k+1} x_{k+1}
  Bi = B(LUdata.p{i},:) - A(LUdata.p{i},LUdata.q{k+1}) * Xhat{k+1};
  % x_i = A_ii \ b_i
  Xhat{i} = LUdata.U{i} \ (LUdata.L{i} \ Bi);
end

% Assemble solution in X.
for i= 1:k
  X(LUdata.q{i},:) = Xhat{i};
end
X(LUdata.q{k+1},:) = Xhat{k+1};
