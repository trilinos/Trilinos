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
  Bi = B(LUdata.Sr{i},:) - A(LUdata.Sr{i},LUdata.c{i}) *(A(LUdata.r{i},LUdata.c{i}) \ B(LUdata.r{i},:));
  Bhat = [Bhat; Bi];
end

% Solve Schur complement.
Xhat{k+1} = LUdata.S{k+1} \ Bhat;

% Solve for diagonal blocks.
for i= 1:k
  % b_i = b_i - A_{i,k+1} x_{k+1}
  Bi = B(LUdata.r{i},:) - A(LUdata.r{i},LUdata.c{k+1}) * Xhat{k+1};
  % x_i = A_ii \ b_i
  Xhat{i} = A(LUdata.r{i},LUdata.c{i}) \ Bi;
end

% Assemble solution in X.
%start=1;
for i= 1:k
  X(LUdata.c{i},:) = Xhat{i};
  %X(LUdata.r{i},:) = Xhat{i};
  %startnext = start+length(LUdata.Sr{i});
  %X(LUdata.Sr{i},:) = Xhat{k+1}(start:startnext-1,:);
  %start = startnext;
end
X(LUdata.c{k+1},:) = Xhat{k+1};
