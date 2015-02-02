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
[partr, partc, permr, permc, domain] = SBBD(A, k);
% Save permutation to SBBD form from partitioning
LUdata.Spermr = permr;
LUdata.Spermc = permc;
% Save the connectivity graph of the partitions
LUdata.domain = domain;

% Save explicit permutation to DBBD form after factoring diagonal blocks and
% Schur complement.
LUdata.Dpermr = [];
LUdata.Dpermc = [];

LUdata.nblocks = k;

% Store original matrix (TODO: we probably only need store parts of A)
LUdata.A = A;

Sc = find(partc==k+1); % columns in local Schur comp.

% Factor diagonal blocks (could be done in parallel)
for i=1:k
  r = find(partr==i); % rows in block i
  c = find(partc==i); % columns in block i
  [nc, dummy] = size(c);
  % Check more rows than columns
  if (length(r) < length(c))
    error('Diagonal block has less rows than columns!');
  end

  % LU factorization on (rectangular) diagonal blocks
  [L, U, p, q] = lu(A(r,c), 'vector');
  pr = r(p(1:nc));     % pivot rows from the rectangular LU
  Sr = r(p(nc+1:end)); % rows in local Schur complement
  LUdata.c{i} = c;
  LUdata.r{i} = pr;
  LUdata.Sr{i} = Sr;
  LUdata.L{i} = L(1:nc, :); % store the square L and throw away the rest.
  LUdata.U{i} = U;  % U is square
  LUdata.p{i} = pr; % same as r{i} when dmperm is not used.
  LUdata.q{i} = c(q);
  LUdata.Dpermr = [LUdata.Dpermr; LUdata.p{i}];
  LUdata.Dpermc = [LUdata.Dpermc; LUdata.q{i}];

  % Form local Schur complement
  % S = S - (A_ki * inv(A_i) * A_ik), where A_i(p,q) = LU
  LUdata.S{i} = A(Sr,Sc) - (A(Sr,c(q))/U) * (L(1:nc, :)\A(pr,Sc)); % Optimized 
end
LUdata.c{k+1} = Sc;

% Assemble and factor global Schur complement
S = [];
Sr = [];
for i=1:k
  S = [S; LUdata.S{i}];
  Sr = [Sr; LUdata.Sr{i}];
end
% Store global S (TODO: Not required if we have all the local S{i}?)
LUdata.S{k+1} = S;

[Snr, dummy] = size(S);
% Assign part numbers to rows of S. Used only in ilu.
LUdata.Spartr = zeros(Snr, 1);
crow = 0;
% TODO: Can extract this from partr
for i=1:k
  [Sinr, dummy] = size(LUdata.S{i});
  LUdata.Spartr(crow+1 : crow+Sinr) = i;
  crow = crow + Sinr;
end

% This is not required, if LUdriver is used. L and U needed if LUsolve is called
% after LUfactor.
i= k+1; % Store Schur complement factors as block k+1
[LUdata.L{i}, LUdata.U{i}, p, q] = lu(S, 'vector');
LUdata.p{i} = p';  % Using local index for the Schur complement rows
LUdata.q{i} = Sc(q);
LUdata.Dpermr = [LUdata.Dpermr; Sr(p)]; % Using global index for the permutation
LUdata.Dpermc = [LUdata.Dpermc; LUdata.q{i}];

% Compute stats
nnzlu = 0;
for i=1:k+1
  nnzlu = nnzlu + nnz(LUdata.L{i}) + nnz(LUdata.U{i});
end
LUdata.nnzlu = nnzlu;
