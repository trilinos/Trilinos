function [part,perm] = zoltan(A, p, method, dir, opt)
% function [part,perm] = zoltan(A, p, method, dir, opt)
%
% Partition a sparse matrix along rows or cols,
% using the Zoltan graph and hypergraph partitioner.
% Balance number of row/columns, or the nonzeros. 
%
% Input: (only first argument is required)
%   A, a sparse matrix
%   p, the number of parts (partitions)
%   method, 'graph' or 'hypergraph'
%   dir, the cut direction (1=rows, 2=cols)
%   opt, Zoltan partitioning options
%
% Output:
%   part, a partition vector with values in 0..p-1
%   perm, a permutation vector to order A by the partition numbers
%
% Written by Erik Boman.
% (C) Copyright Sandia Corporation, 2006-2007

if (nargin<5)
  opt = [];
end
if (nargin<4)
  dir = 1; 
end
if (nargin<3)
  method = 'hypergraph';
end
if (nargin<2)
  p = 2;
end

% For graph model, symmetrize
if (strcmp(method,'graph'))
  [m,n] = size(A);
  if  (m==n) % Square
    S = (A~=0); % Structure only
    if (norm(S-S','fro'))
      A = A+A';
    end
  else % Rectangular 
    if (dir==1)
      A = A*A';
    else
      A = A'*A;
    end
  end
end

% Write matrix to file
mmwrite('matlab.mtx', A);

% Copy standard zdrive input file  (overwrite old zdrive.inp)
copyfile ('zdrive.matlab', 'zdrive.inp');
fp = fopen('zdrive.inp', 'a');

% Append load balance method
fprintf(fp, 'Decomposition method    = %s\n', method);

%% Direction determines row or column partition
if (dir==1)
  fprintf(fp, 'File Type               = matrixmarket, objects=rows\n');
elseif (dir==2)
  fprintf(fp, 'File Type               = matrixmarket, objects=cols\n');
else
  error('Invalid value for dir; must be 1 or 2!');
end

% Append number of parts
fprintf(fp, 'Zoltan parameter = num_global_partitions=%d\n', p);
% Append other options
if (opt)
  % Loop over options
  for option = opt
    fprintf(fp, 'Zoltan parameter = %s\n', option);
  end
end

% Run zdrive to partition the matrix
% zdrive must be in your path.
!zdrive

% Parse zdrive output file  
fp = fopen('matlab.out.1.0', 'r');
if (~fp)
  error('Could not open zdrive output file\n');
end
% Skip all lines before 'GID'
word = 'abc';
while (~strncmp(word, 'GID', 3))
  % Read only first word, ignore rest of line
  word = fscanf(fp,'%s',1);
  fscanf(fp,'%*[^\n]',1);
end
% Read the partition numbers; file has 4 fields per line
[P, num] = fscanf(fp, '%d', [4,inf]);
% First two rows in P (columns in output file) give the partition vector
part = zeros(size(P(1,:)));
part(P(1,:)) = P(2,:);

% Construct corresponding permutation vector 
perm = [];
for i= 0:p-1
  perm = [perm, find(part==i)];
end

