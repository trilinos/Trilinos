#! /usr/bin/env python
import os
import sys
import re
import Numeric

import setpath
import Epetra
import EpetraExt
import AztecOO
import ML
import IFPACK
import Amesos
import Galeri

# -------------------------------------------------------------------------
def generator(problemID, comm):
  if (problemID == "Laplace2D_100_100"):
    GaleriList = {
      "nx": 100,
      "ny": 100,
      "mx": 1,
      "my": 1
    }

    Map = Galeri.CreateMap("Cartesian2D", comm, GaleriList);

    Matrix = Galeri.CreateCrsMatrix("Laplace2D", Map, GaleriList);

    LHS = Epetra.Vector(Map);
    RHS = Epetra.Vector(Map);
    ExactSolution = Epetra.Vector(Map); ExactSolution.Random();
    Matrix.Apply(ExactSolution, RHS);
    LHS.PutScalar(0.0);

  elif (problemID == "bcsstk13"):
    FileName = "/people_old/trilinos_www/matrices/bcsstk13.rsa";
    Map, Matrix, LHS, RHS, ExactSolution = Galeri.ReadHB(FileName, comm);

  return(Map, Matrix, LHS, RHS, ExactSolution);

# -------------------------------------------------------------------------
def analyze(Map, Matrix, LHS, RHS, ExactSolution):
  IFPACK.AnalyzeMatrix(Matrix);
  IFPACK.AnalyzeMatrixElements(Matrix);

# -------------------------------------------------------------------------
def iterative(Map, Matrix, LHS, RHS, ExactSolution, List):
  Prec = ML.MultiLevelPreconditioner(Matrix, False);
  Prec.SetParameterList(List);
  Prec.ComputePreconditioner();

  Solver = AztecOO.AztecOO(Matrix, LHS, RHS)
  Solver.SetPrecOperator(Prec)
  if (List['solver'] == "gmres"):
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_gmres);
  elif List['solver'] == "cg":
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_gmres);
  else:
    print "Solver type not correct"
  Solver.SetAztecOption(AztecOO.AZ_output, 16);
  err = Solver.Iterate(List['iters'], List['tol']) 

# -------------------------------------------------------------------------
def direct(Map, Matrix, LHS, RHS, ExactSolution, List):
  List['PrintStatus'] = True;
  List['PrintTiming'] = True;

  Problem = Epetra.LinearProblem();
  Problem.SetOperator(Matrix)
  Problem.SetLHS(LHS)
  Problem.SetRHS(RHS);
  Factory = Amesos.Factory();
  Solver = Factory.Create(List['solver'], Problem);
  Solver.SetParameters(List)
  Solver.SymbolicFactorization()
  Solver.NumericFactorization()
  Solver.Solve()

# -------------------------------------------------------------------------
def main():

  # Grab config file from command line ========================================
  
  if sys.argv[1]:
    configFile = sys.argv[1]
  else:
    print "Please supply a config file." ; return    
  
  # Parse config file =========================================================
  
  file = open(configFile, 'r')
  fileContents = file.read()
  pattern = re.compile(r"^(\w+)\s*=(.*)$", re.M)
  config = pattern.findall(fileContents)
  config = dict(config)
  for key in config:
    config[key] = config[key].strip()
  
  # Parse list ================================================================
  
  List = {}

  problemID = config['PROBLEM_ID'];
  step = config['STEP'];

  if step == 'step_2':
    counter = 0;
    donothing = 0;
  elif step == 'step_3i':
    counter = int(config['COUNTER']);
    List['solver'] = config['SOLVER'];
    List['iters'] = int(config['ITERS']);
    List['tol'] = float(config['TOL']);
  elif step == 'step_3d':
    counter = int(config['COUNTER']);
    List['solver'] = config['SOLVER'];

  for i in range(1,counter-1):
    name = config['__PyTrilinos__name' + str(i)];
    type = config['__PyTrilinos__type' + str(i)];
    value = config['__PyTrilinos__value' + str(i)];

    if (type == 'int'):
      List[name] = int(value);
    elif (type == 'string'):
      List[name] = value;
    elif (type == 'double'):
      List[name] = float(value);
    elif (type == ''):
      continue;
    else:
      print "Type ", type, " not recognized"

  # Construct the problem =====================================================

  comm = Epetra.PyComm()
  
  (Map, Matrix, LHS, RHS, ExactSolution) = generator(problemID, comm);


  if step == 'step_2':
    analyze(Map, Matrix, LHS, RHS, ExactSolution);
  elif step == 'step_3i':
    iterative(Map, Matrix, LHS, RHS, ExactSolution, List);
  elif step == 'step_3d':
    direct(Map, Matrix, LHS, RHS, ExactSolution, List);

# -------------------------------------------------------------------------
if __name__ == "__main__":
  main()
