from PyTrilinos import Epetra, ML, EpetraExt, IFPACK, AztecOO
from math import sin
comm = Epetra.PyComm()
(ierr, Matrix) = EpetraExt.MatrixMarketFileToCrsMatrix("/tmp/MPSalsa_Diamond_med.mm", comm);

Map = Matrix.RowMatrixRowMap()


# -------------------------------------- #
# building solutions and right-hand side #
# -------------------------------------- #

LHS = Epetra.Vector(Map)
RHS = Epetra.Vector(Map)
LHSView = LHS.ExtractView()
n = LHS.MyLength()
for i in xrange(0, n):
  LHSView[i] = sin(i * 3.1415 / n) * sin(i * 3.1415 / n)

Matrix.Multiply(False, LHS, RHS)


# ------------------------------------------------- #
# Parameters to run ML. Remember the RCM reordering #
# for Charon matrices!                              #
# ------------------------------------------------- #

MLList = {
  "max levels"           : 10,
  "output"               : 10,
  "smoother: type"       : "IFPACK",
  "smoother: ifpack type" : "ILU",
  "smoother: ifpack list" : {
    "fact: relax value": 0.0, 
    "fact: absolute threshold": 0e-2,
    "fact: level-of-fill": 0,
    "schwarz: reordering type": "none"
  },
  "smoother: ifpack level-of-fill" : 0,
  "smoother: pre or post": "post",
  "aggregation: use tentative restriction": False,
  "PDE equations"        : 3,
  "aggregation: damping factor" : 0.0,
  "aggregation: threshold": 0.00,
  "aggregation: type"    : "Uncoupled",
  "aggregation: nodes per aggregate": 64,
  "aggregation: block scaling": False,
  "energy minimization: enable": False,
  "energy minimization: type": 2
}

Prec = ML.MultiLevelPreconditioner(Matrix, False)
Prec.SetParameterList(MLList)
Prec.ComputePreconditioner()

Solver = AztecOO.AztecOO(Matrix, LHS, RHS)
Solver.SetPrecOperator(Prec)
Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_gmres);
Solver.SetAztecOption(AztecOO.AZ_kspace, 200)
Solver.SetAztecOption(AztecOO.AZ_output, 16);
LHS.PutScalar(0.0)
err = Solver.Iterate(155, 1e-4)


# ----------------------------------------------- #
# parameters to run one-level schemes with ifpack #
# ----------------------------------------------- #

Factory = IFPACK.Factory()
Prec = Factory.Create("ILU", Matrix);
IFPACKList = {
  "fact: absolute threshold": 0e-9,
  "fact: relative threshold": 1.0,
  "fact: relax value": 0.0, 
  "fact: level-of-fill": 0,
  "schwarz: reordering type": "none"
  }
Prec.SetParameters(IFPACKList)
Prec.Initialize() 
Prec.Compute()
LHS.PutScalar(0.0)
Solver = AztecOO.AztecOO(Matrix, LHS, RHS)
Solver.SetAztecOption(AztecOO.AZ_kspace, 200)
Solver.SetPrecOperator(Prec)
Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_gmres)
Solver.SetAztecOption(AztecOO.AZ_output, 16)
Solver.Iterate(255, 1e-4)




