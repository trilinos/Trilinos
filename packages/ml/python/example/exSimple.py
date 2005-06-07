#! /usr/bin/env python
try:
  import setpath
  import ML, Triutils, AztecOO, Epetra
except:
  try:
    from PyTrilinos import ML, Triutils, AztecOO, Epetra
  except ImportError:
    raise ImportError, "error w/ ML or Triutils or AztecOO or Epetra"

# builds the linear system matrix and sets up starting solution and
# right-hand side
nx = 100;
ny = 100;
Comm = Epetra.SerialComm();
Gallery = Triutils.CrsMatrixGallery("laplace_2d", Comm)
Gallery.Set("nx", nx);
Gallery.Set("ny", ny);
Matrix = Gallery.GetMatrix();
LHS = Gallery.GetStartingSolution();
RHS = Gallery.GetRHS();

# sets up the parameters for ML using a python dictionary
MLList = {
  "max levels"        : ("int", "3"), 
  "output"            : ("int", "10"),
  "smoother: type"    : ("string", "symmetric Gauss-Seidel"),
  "aggregation: type" : ("string", "Uncoupled")
};

# creates the preconditioner and computes it
Prec = ML.MultiLevelPreconditioner(Matrix, False);
Prec.SetParameterList(MLList);
Prec.ComputePreconditioner();

# sets up the solver, specifies Prec as preconditioner, and
# solves using CG.
Solver = AztecOO.AztecOO(Matrix, LHS, RHS)
Solver.SetPrecOperator(Prec)
Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg);
Solver.SetAztecOption(AztecOO.AZ_output, 16);
Solver.Iterate(1550, 1e-5)
