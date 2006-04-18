#! /usr/bin/env python
try:
    import setpath
except ImportError:
    from PyTrilinos import Epetra, AztecOO, Triutils, ML
    print "Using installed versions of ML, Triutils, AztecOO, Epetra"
else:
    import Epetra, AztecOO, Triutils, ML

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
