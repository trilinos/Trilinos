
# this script is supposed to be cut-and-paste'd into python, as a way to
# interactively test matrices. 

from PyTrilinos import  Epetra, AztecOO, ML, Galeri
from Numeric import sin

import sys

def runit(size, prec, diff):
  # builds the linear system matrix and sets up starting solution and
  # right-hand side
  Comm = Epetra.PyComm()

  print "SIZE = ", size
  print "PREC = ", prec
  print "DIFF = ", diff
  
  GaleriList = {
    "nx": size,
    "ny": size, 
    "mx": 1, 
    "my": 1,
    "conv": 1.0,
    "diff": diff
  }
  
  Map = Galeri.CreateMap("Cartesian2D", Comm, GaleriList)
  Matrix = Galeri.CreateCrsMatrix("BentPipe2D", Map, GaleriList)
  
  Util = Epetra.Util()
  Util.SetSeed(0)
  LHS = Epetra.Vector(Map)
  RHS = Epetra.Vector(Map)
  n = LHS.MyLength()
  for i in xrange(0, n):
    LHS[i] = sin(i * 3.1415 / n) * sin(i * 3.1415 / n)
  Matrix.Multiply(False, LHS, RHS)
  LHS.PutScalar(0.0)
  
  # sets up the parameters for ML using a python dictionary
  if prec == "SA":
    MLList = {
      "max levels"        : 10, 
      "output"            : 10,
      "increasing or decreasing": "increasing",
      "aggregation: type" : "Uncoupled-MIS",
      "smoother: type (level 0)"    : "symmetric Gauss-Seidel",
      "smoother: damping factor": 0.67,
      "smoother: sweeps": 1,
      "smoother: pre or post": "both",
      "PDE equations": 1,
      "aggregation: damping factor": 1.333,
      "eigen-analysis: type": "power-method"
    }
  elif prec == "NSA":
    MLList = {
      "max levels"        : 10, 
      "output"            : 10,
      "increasing or decreasing": "increasing",
      "aggregation: type" : "Uncoupled-MIS",
      "smoother: type (level 0)"    : "symmetric Gauss-Seidel",
      "smoother: damping factor": 0.67,
      "smoother: sweeps": 1,
      "smoother: pre or post": "both",
      "PDE equations": 1,
      "aggregation: damping factor": 0.0000,
      "eigen-analysis: type": "power-method"
    }
  elif prec == "NSR":
    MLList = {
      "max levels"        : 10, 
      "output"            : 10,
      "increasing or decreasing": "increasing",
      "aggregation: type" : "Uncoupled-MIS",
      "smoother: type (level 0)"    : "symmetric Gauss-Seidel",
      "smoother: damping factor": 0.67,
      "smoother: sweeps": 1,
      "smoother: pre or post": "both",
      "PDE equations": 1,
      "aggregation: use tentative restriction": True,
      "aggregation: damping factor": 1.0000,
      "eigen-analysis: type": "power-method"
    }
  elif prec == "EMIN":
    MLList = {
      "max levels"        : 10, 
      "output"            : 10,
      "increasing or decreasing": "increasing",
      "aggregation: type" : "Uncoupled-MIS",
      "smoother: type (level 0)"    : "symmetric Gauss-Seidel",
      "smoother: damping factor": 0.67,
      "smoother: sweeps": 1,
      "smoother: pre or post": "both",
      "PDE equations": 1,
      "energy minimization: enable": True,
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
  Solver.SetAztecOption(AztecOO.AZ_output, 32);
  Solver.Iterate(400, 1e-8)
  
  del Prec

if __name__ == "__main__":

  #runit(64,   "NSA", 1e-5);
  #runit(128,  "EMIN", 1e-5);
  #runit(256,  "EMIN", 1e-5);
  runit(512,  "SA", 1e-5);
  runit(1024, "SA", 1e-5);
