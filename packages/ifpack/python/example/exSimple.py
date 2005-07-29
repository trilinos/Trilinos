#! /usr/bin/env python
try:
  import setpath
  import Epetra, Triutils, AztecOO, IFPACK
except:
  from PyTrilinos import Epetra, Triutils, AztecOO, IFPACK
  print "Using system-installed Epetra, Triutils, AztecOO, IFPACK"

import sys

def main():

  # read the matrix from file in H/B format. The filename is specified by the
  # first argument in the compile line. If no filename is specified, then the
  # code build a matrix using matrix gallery.
  Comm = Epetra.PyComm();

  args = sys.argv[1:]
  if len(args) == 0:
    Gallery = Triutils.CrsMatrixGallery("recirc_2d", Comm)
    Gallery.Set("nx", 100)
    Gallery.Set("ny", 100)
    Matrix = Gallery.GetMatrix()
    Exact = Gallery.GetExactSolution()
    LHS = Gallery.GetStartingSolution()
    RHS = Gallery.GetRHS()
  else:
    Map, Matrix, LHS, RHS, Exact = Triutils.ReadHB(args[0], Comm);

  # Creates the IFPACK preconditioner, in this case an incomplete
  # Cholesky factorization, with fill-in of 5
  Factory = IFPACK.Factory();
  Prec = Factory.Create("ILU", Matrix);
  IFPACKList = {
    "fact: level-of-fill": 1
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


if __name__ == "__main__":
  main()
