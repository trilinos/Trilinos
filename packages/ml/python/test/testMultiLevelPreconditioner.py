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
    import AztecOO
    import Triutils
    import ML
else:
    try:
        import setpath
        import Epetra
        import AztecOO
        import Triutils
        import ML
    except ImportError:
        from PyTrilinos import Epetra
        from PyTrilinos import AztecOO
        from PyTrilinos import Triutils
        from PyTrilinos import ML
        print >>sys.stderr, "Using installed versions of ML, Triutils, AztecOO, Epetra"

# builds the linear system matrix and sets up starting solution and
# right-hand side
nx = 100
ny = 100
Comm = Epetra.PyComm()
Gallery = Triutils.CrsMatrixGallery("laplace_2d", Comm)
Gallery.Set("nx", nx)
Gallery.Set("ny", ny)
Matrix = Gallery.GetMatrix()
LHS = Gallery.GetStartingSolution()
RHS = Gallery.GetRHS()

# sets up the parameters for ML using a python dictionary
MLList = {"max levels"        : 3, 
          "output"            : 10,
          "smoother: type"    : "symmetric Gauss-Seidel",
          "aggregation: type" : "Uncoupled"             }

# creates the preconditioner and computes it
Prec = ML.MultiLevelPreconditioner(Matrix, False)
Prec.SetParameterList(MLList)
Prec.ComputePreconditioner()

# sets up the solver, specifies Prec as preconditioner, and
# solves using CG.
Solver = AztecOO.AztecOO(Matrix, LHS, RHS)
Solver.SetPrecOperator(Prec)
Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg)
Solver.SetAztecOption(AztecOO.AZ_output, 16)
Solver.Iterate(1550, 1e-5)

if Comm.MyPID() == 0: print "End Result: TEST PASSED"
