function LUdata = LUfactor(A, k)
% function LUdata = LUfactor(A, k)
%
% 2-level sparse LU based on SBBD form (from hypergraph 
% partitioning). Factorization phase only. Use LUsolve to solve.
%
% Input:
%   A, matrix to be factored.
%   k, (integer) number of parts.
%
% Erik Boman, Sandia National Labs, 2010.

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

% Singly bordered block diagonal form.
[partr, partc] = SBBD(A, k);

LUdata.nblocks = k;

% Store original matrix (TODO: we probably only need store parts of A)
LUdata.A = A;

% Factor diagonal blocks (could be done in parallel)
for i=1:k
  r = find(partr==i); % rows in block i
  c = find(partc==i); % columns in block i
  LUdata.c{i} = c;
  % Check more rows than columns
  if (length(r) < length(c))
    error('Diagonal block has less rows than columns!');
  end
  % dmperm matching to get zero-free diagonal (avoid rank-deficiency)
  dm = dmperm(A(r,c));
  Sr = r; Sr(dm) = []; % rows in local Schur complement
  LUdata.Sr{i} = Sr;
  r = r(dm);           % rows to eliminate
  LUdata.r{i} = r;
  % LU factorization on (square) diagonal blocks
  [L, U, p, q] = lu(A(r,c), 'vector');
  LUdata.L{i} = L;
  LUdata.U{i} = U;
  LUdata.p{i} = p;
  LUdata.q{i} = q;
  % Form local Schur complement
  Sc = find(partc==k+1); % columns in local Schur comp.
  % S = S - (A_ki * inv(A_i) * A_ik), where A_i(p,q) = LU
  LUdata.S{i} = A(Sr,Sc) - A(Sr,c) * (A(r,c) \ A(r,Sc)); % Naive method
  %LUdata.S{i} = A(Sr,Sc) - (A(Sr,c(q))/U) * (L\A(r(p),Sc)); % Optimized 
end
LUdata.c{k+1} = Sc;

% Assemble and factor global Schur complement
S = [];
for i=1:k
  S = [S; LUdata.S{i}];
end
% Store global S (TODO: Not required if we have all the local S{i}?)
LUdata.S{k+1} = S;

i= k+1; % Store Schur complement factors as block k+1
[LUdata.L{i}, LUdata.U{i}, LUdata.p{i}, LUdata.q{i}] = lu(S);

% Compute stats
nnzlu = 0;
for i=1:k+1
  nnzlu = nnzlu + nnz(LUdata.L{i}) + nnz(LUdata.U{i});
end
LUdata.nnzlu = nnzlu;
