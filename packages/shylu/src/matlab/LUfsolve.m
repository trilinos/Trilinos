function Bhat = LUfsolve(LUdata, B)
% function Bhat = LUfsolve(LUdata, B)
%
% Compute the rhs corresponding to the Schur complement.
% LUfactor must be called first. Schur complement must be solved and then
% LUbsolve should be used.
%
% Erik Boman, Sandia National Labs, 2010.
% Siva Rajamanickam, Sandia National Labs, 2010.

k = LUdata.nblocks;
A = LUdata.A;

% Modify rhs for Schur complement.
Bhat = [];
for i= 1:k
  % b_i = b_i - \bar{A_i} * (A_i \ b_i)
  Bi = B(LUdata.Sr{i},:) - (A(LUdata.Sr{i},LUdata.q{i})/LUdata.U{i}) * ...
		(LUdata.L{i}\ B(LUdata.r{i},:));
  Bhat = [Bhat; Bi];
end
