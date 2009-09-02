# this script is supposed to be cut-and-paste'd into python, as a way to
# interactively test matrices. 

from PyTrilinos import  Epetra, EpetraExt, AztecOO, ML

from Numeric import sin

# builds the linear system matrix and sets up starting solution and
# right-hand side
Comm = Epetra.PyComm()

# I got similar results for _1 and _4
MapName = "/Users/marzio/matrices/PREMO/Falcon-step2/Falcon_ss_1_map"
MatrixName = "/Users/marzio/matrices/PREMO/Falcon-step2/Falcon_ss_4_matrix"
LHSName = "/Users/marzio/matrices/PREMO/Falcon-step2/Falcon_ss_1_guess"
RHSName = "/Users/marzio/matrices/PREMO/Falcon-step2/Falcon_ss_1_guess"

(ierr, Map) = EpetraExt.MatrixMarketFileToBlockMap(MapName, Comm)
(ierr, LHS) = EpetraExt.MatrixMarketFileToMultiVector(LHSName, Map)
(ierr, RHS) = EpetraExt.MatrixMarketFileToMultiVector(RHSName, Map)
(ierr, Matrix) = EpetraExt.MatrixMarketFileToCrsMatrix(MatrixName, Map)

LHSView = LHS.ExtractView()
n = LHS.MyLength()
for i in xrange(0, n):
  LHSView[0][i] = sin(i * 3.1415 / n) * sin(i * 3.1415 / n)

Matrix.Multiply(False, LHS, RHS)
LHS.PutScalar(0.0)

# sets up the parameters for ML using a python dictionary
MLList = {
  "max levels"        : 10, 
  "output"            : 10,
  "increasing or decreasing": "increasing",
  "aggregation: type" : "Uncoupled",
  "smoother: type (level 0)"    : "Aztec",
  "smoother: damping factor": 0.67,
  "smoother: sweeps": 1,
  "smoother: pre or post": "post",
  "smoother: type (level 1)"    : "Aztec",
  "PDE equations": 5,
  "aggregation: use tentative restriction": False,
  "aggregation: damping factor": 0.000,
  "eigen-analysis: type": "power-method",
  "aggregation: block scaling": False,
  "energy minimization: enable": False,
  "energy minimization: type": 2
}

# creates the preconditioner and computes it
Prec = ML.MultiLevelPreconditioner(Matrix, False)
Prec.SetParameterList(MLList)
Prec.ComputePreconditioner()

# sets up the solver, specifies Prec as preconditioner, and
# solves using CG.
Solver = AztecOO.AztecOO(Matrix, LHS, RHS)
Solver.SetPrecOperator(Prec)
Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_gmres);
Solver.SetAztecOption(AztecOO.AZ_kspace, 30);
Solver.SetAztecOption(AztecOO.AZ_output, 1);
Solver.Iterate(200, 1e-8)

del Prec
