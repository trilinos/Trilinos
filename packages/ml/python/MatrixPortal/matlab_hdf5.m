function write_hdf5(FileName, A, RHS, LHS, ExactSolution, NullSpace)
% Usage: WRITE_HDF5(FileName, A, RHS, LHS, ExactSolution, NullSpace)
%
% - A: the linear system matrix. Must be square and sparse
%
% - RHS: the right-hand side vector. 
%
% - LHS: the starting solution vector.
%
% - ExactSolution: the exact solution of the linear system
%
% - NullSpace: an approximation of the null space, to be used
%              by smoothed aggregation preconditioners.
%
% Author: Marzio Sala, ETHZ/COLAB
%
% Date: Last updated on 04-May-06

setLHS = 1;
setRHS = 1;
setExactSolution = 1;
setNullSpace = 1;

if (nargin == 2)
  setLHS = 0;
  setRHS = 0;
  setExactSolution = 0;
  setNullSpace = 0;
elseif (nargin == 3)
  setRHS = 0;
  setExactSolution = 0;
  setNullSpace = 0;
elseif (nargin == 4)
  setExactSolution = 0;
  setNullSpace = 0;
elseif (nargin == 5)
  setNullSpace = 0;
end

n = size(A, 1);
if n ~= size(A, 2)
  disp('Matrix must be square');
  return
end

[ROW,COL,VAL] = find(A);
hdf5write(FileName, '/MATRIX/__type__',           'Epetra_RowMatrix');
hdf5write(FileName, '/MATRIX/NumGlobalRows',      int32(n),       'WriteMode', 'append');
hdf5write(FileName, '/MATRIX/NumGlobalCols',      int32(n),       'WriteMode', 'append');
hdf5write(FileName, '/MATRIX/NumGlobalNonzeros',  int32(n),       'WriteMode', 'append');
hdf5write(FileName, '/MATRIX/NumGlobalDiagonals', int32(n),       'WriteMode', 'append');
hdf5write(FileName, '/MATRIX/MaxNumEntries',      int32(1),       'WriteMode', 'append');
hdf5write(FileName, '/MATRIX/NormOne',            1.0,            'WriteMode', 'append');
hdf5write(FileName, '/MATRIX/NormInf',            1.0,            'WriteMode', 'append');
hdf5write(FileName, '/MATRIX/ROW',                int32(ROW - 1), 'WriteMode', 'append');
hdf5write(FileName, '/MATRIX/COL',                int32(COL - 1), 'WriteMode', 'append');
hdf5write(FileName, '/MATRIX/VAL',                VAL,            'WriteMode', 'append');

if setRHS == 1
  if n ~= size(RHS, 1)
    disp('RHS has a different dimension than A');
    disp('All vectors must be row vectors');
    return
  end

  if size(RHS, 2) > 1
    disp('RHS is a multivector; only one RHS is supported');
    return
  end

  hdf5write(FileName, '/RHS/__type__',    'Epetra_MultiVector', 'WriteMode', 'append');
  hdf5write(FileName, '/RHS/GlobalLength',int32(n),             'WriteMode', 'append');
  hdf5write(FileName, '/RHS/NumVectors',  int32(1),             'WriteMode', 'append');
  hdf5write(FileName, '/RHS/Values',      RHS,                  'WriteMode', 'append');
end

if setLHS == 1
  if n ~= size(LHS, 1)
    disp('LHS has a different dimension than A');
    disp('All vectors must be row vectors');
    return
  end

  if size(LHS, 2) > 1
      disp('LHS is a multivector; only one LHS is supported');
    return
  end

  hdf5write(FileName, '/LHS/__type__',    'Epetra_MultiVector', 'WriteMode', 'append');
  hdf5write(FileName, '/LHS/GlobalLength',int32(n),             'WriteMode', 'append');
  hdf5write(FileName, '/LHS/NumVectors',  int32(1),             'WriteMode', 'append');
  hdf5write(FileName, '/LHS/Values',      LHS,                  'WriteMode', 'append');
end

if setExactSolution == 1
  if n ~= size(ExactSolution, 1)
    disp('ExactSolution has a different dimension than A');
    disp('All vectors must be row vectors');
    return
  end

  if size(ExactSolution, 2) > 1
      disp('ExactSolution is a multivector; only one exact solution  is supported');
    return
  end

  hdf5write(FileName, '/ExactSolution/__type__',    'Epetra_MultiVector', 'WriteMode', 'append');
  hdf5write(FileName, '/ExactSolution/GlobalLength',int32(n),             'WriteMode', 'append');
  hdf5write(FileName, '/ExactSolution/NumVectors',  int32(1),             'WriteMode', 'append');
  hdf5write(FileName, '/ExactSolution/Values',      ExactSolution,        'WriteMode', 'append');
end

if setNullSpace == 1
  if n ~= size(NullSpace, 1)
    disp('ExactSolution has a different dimension than A');
    disp('All vectors must be row vectors');
    return
  end

  NumVectors = size(NullSpace, 2);
  hdf5write(FileName, '/NULLSPACE/__type__',    'Epetra_MultiVector', 'WriteMode', 'append');
  hdf5write(FileName, '/NULLSPACE/GlobalLength',int32(n),             'WriteMode', 'append');
  hdf5write(FileName, '/NULLSPACE/NumVectors',  int32(NumVectors),    'WriteMode', 'append');
  hdf5write(FileName, '/NULLSPACE/Values',      ExactSolution,        'WriteMode', 'append');
end
