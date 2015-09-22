function X = LUbsolve(LUdata, B)
% function X = LUbsolve(LUdata, B)
%
% Solve for the diagonal blocks from LUfactor.
% LUfactor must be called first. Schur complement must be solved and then
% LUbsolve should be used.
%
% Erik Boman, Sandia National Labs, 2010.
% Siva Rajamanickam, Sandia National Labs, 2010.

k = LUdata.nblocks;
A = LUdata.A;
Xhat{k+1} = LUdata.Xhat; % Solution for the Schur complement

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
