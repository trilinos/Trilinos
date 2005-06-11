#! /usr/bin/env python
try:
  from PyTrilinos import IFPACK, Triutils, AztecOO, Epetra
except:
  raise ImportError, "error w/ IFPACK or Triutils or AztecOO or Epetra"


# read the matrix from file, here `bcsstk01.rsa' in HB format
Epetra.Init()
Comm = Epetra.PyComm();
Map, Matrix, LHS, RHS, Exact = Triutils.ReadHB("bcsstk01.rsa", Comm);

# Creates the IFPACK preconditioner, in this case an incomplete
# Cholesky factorization, with fill-in of 5
Factory = IFPACK.Factory();
Prec = Factory.Create("IC", Matrix);
IFPACKList = {
  "fact: level-of-fill": ("int", "5")
}
Prec.SetParameters(IFPACKList);
Prec.Initialize();
Prec.Compute();

# Creates the AztecOO solver, using GMRES and IFPACK as preconditioner
Problem = Epetra.LinearProblem(Matrix, LHS, RHS);
Solver = AztecOO.AztecOO(Problem);
Solver.SetPrecOperator(Prec)
Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_gmres);
Solver.SetAztecOption(AztecOO.AZ_output, 16);
Solver.Iterate(1550, 1e-5)

print Prec

# Computes the 2-norm of the true residual
LHS.Update(1.0, Exact, -1.0);
norm = LHS.Norm2()[1] / RHS.Norm2()[1]
if Comm.MyPID() == 0:
  print "After solution of the linear system:"
  print "||x - x_exact||_2 / ||b||_2 = ", norm

Epetra.Finalize()
