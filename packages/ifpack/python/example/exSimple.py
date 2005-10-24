#! /usr/bin/env python
import sys

try:
  import setpath
  import Epetra
  import Galeri
  import AztecOO
  import IFPACK
except:
  try: 
    # Loading Epetra should ensure MPI_Init() and MPI_Finalize()
    from PyTrilinos import Epetra
    from PyTrilinos import Galeri, AztecOO, IFPACK
    print "Using system-installed Epetra, Galeri, AztecOO, IFPACK"
  except:
    print "You need the Epetra, AztecOO, Galeri and IFPACK"
    print "modules to run this examples."
    sys.exit(0)

def main():

  # read the matrix from file in H/B format. The filename is specified by the
  # first argument in the compile line. If no filename is specified, then the
  # code build a matrix using matrix gallery.
  Comm = Epetra.PyComm();

  args = sys.argv[1:]
  if len(args) == 0:
    nx = 30; ny = 30
    GaleriList = {
      "n": nx * ny,
      "nx": nx,
      "ny": ny
    }
    Map = Galeri.CreateMap("Linear", Comm, GaleriList)
    Matrix = Galeri.CreateCrsMatrix("Recirc2D", Map, GaleriList)
    Exact = Epetra.Vector(Map); 
    LHS = Epetra.Vector(Map); 
    RHS = Epetra.Vector(Map); 
    Exact.Random()       # fix exact solution
    LHS.PutScalar(0.0)   # fix starting solution
    Matrix.Multiply(False, Exact, RHS) # fix rhs corresponding to Exact
  else:
    try:
      Map, Matrix, LHS, RHS, Exact = Triutils.ReadHB(args[0], Comm);
    except:
      print "Specified matrix cannot be found"
      exit(0)

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
