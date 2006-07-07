#! /usr/bin/env python

from   optparse import *
import sys

parser = OptionParser()
parser.add_option("-t", "--testharness", action="store_true",
                  dest="testharness", default=False,
                  help="test local build modules; prevent loading system-installed modules")
parser.add_option("-v", "--verbosity", type="int", dest="verbosity", default=2,
                  help="set the verbosity level [default 2]")
options,args = parser.parse_args()
if options.testharness:
  import setpath
  import Epetra
  import Triutils
  import Galeri
  import AztecOO
  import IFPACK
else:
  try:
    import setpath
    import Epetra
    import Triutils
    import Galeri
    import AztecOO
    import IFPACK
  except:
    from PyTrilinos import Epetra
    from PyTrilinos import Triutils
    from PyTrilinos import Galeri
    from PyTrilinos import AztecOO
    from PyTrilinos import IFPACK
    print >>sys.stderr, "Using system-installed Epetra, Galeri, AztecOO, IFPACK"

comm    = Epetra.PyComm()
iAmRoot = comm.MyPID() == 0

def main():

  # read the matrix from file in H/B format. The filename is specified by the
  # first argument in the compile line. If no filename is specified, then the
  # code build a matrix using matrix gallery.
  if len(args) == 0:
    nx = 30; ny = 30
    GaleriList = {"n"  : nx * ny,
                  "nx" : nx,
                  "ny" : ny
                  }
    Map = Galeri.CreateMap("Linear", comm, GaleriList)
    Matrix = Galeri.CreateCrsMatrix("Recirc2D", Map, GaleriList)
    Exact = Epetra.Vector(Map) 
    LHS = Epetra.Vector(Map) 
    RHS = Epetra.Vector(Map) 
    Exact.Random()       # fix exact solution
    LHS.PutScalar(0.0)   # fix starting solution
    Matrix.Multiply(False, Exact, RHS) # fix rhs corresponding to Exact
  else:
    Map, Matrix, LHS, RHS, Exact = Triutils.ReadHB(args[0], comm)

  # Creates the IFPACK preconditioner, in this case an incomplete
  # Cholesky factorization, with fill-in of 5
  Factory = IFPACK.Factory()
  Prec = Factory.Create("ILU", Matrix)
  IFPACKList = {
    "fact: level-of-fill": 1
    }
  Prec.SetParameters(IFPACKList)
  Prec.Initialize()
  Prec.Compute()

  # Creates the AztecOO solver, using GMRES and IFPACK as preconditioner
  Problem = Epetra.LinearProblem(Matrix, LHS, RHS)
  Solver = AztecOO.AztecOO(Problem)
  Solver.SetPrecOperator(Prec)
  Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_gmres)
  Solver.SetAztecOption(AztecOO.AZ_output, 16)
  Solver.Iterate(1550, 1e-5)

  print Prec

  # Computes the 2-norm of the true residual
  LHS.Update(1.0, Exact, -1.0)
  norm = LHS.Norm2() / RHS.Norm2()
  if iAmRoot:
    print "After solution of the linear system:"
    print "||x - x_exact||_2 / ||b||_2 = ", norm

if __name__ == "__main__":
  main()
  if iAmRoot: print "End Result: TEST PASSED"
