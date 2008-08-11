#! /usr/bin/env python
import os
import sys
import re
import string
import Numeric

import setpath
import Epetra
import EpetraExt
import AztecOO
import ML
import IFPACK
import Amesos
import Galeri

comm = Epetra.PyComm()
  
def set_type(List, name, type, value):  
  if (type == 'int'):
    List[name] = int(value);
  elif (type == 'string'):
    List[name] = value;
  elif (type == 'double'):
    List[name] = float(value);
  elif (type == 'bool'):
    List[name] = bool(value);
  elif (type == ''):
    donothing = ""
  else:
    print "Type ", type, " not recognized"

# -------------------------------------------------------------------------
def generator(problemID, comm):
  GaleriList = {}
  if problemID[0:3] == "MM_":
    FileName = "/people_old/trilinos_www/matrices/" + problemID[3:];
    Map, Matrix, LHS, RHS, ExactSolution = Galeri.ReadHB(FileName, comm);

  else:
    parts = string.split(problemID, '_');
    ProblemType = parts[0];
    for i in range(1, len(parts)):
      p2 = string.split(parts[i], '=')
      type = p2[0][0];
      name = p2[0][1:]
      value = p2[1];
      if (type == "i"):
        GaleriList[name] = int(value);
      elif type == "d":
        GaleriList[name] = float(value);
      elif type == "s":
        GaleriList[name] = value;
    
    if string.find(ProblemType, '2D') != -1:
      MapType = "Cartesian2D";
    elif string.find(ProblemType, '3D') != -1:
      MapType = "Cartesian3D"

    Map = Galeri.CreateMap(MapType, comm, GaleriList);

    Matrix = Galeri.CreateCrsMatrix(ProblemType, Map, GaleriList);

    LHS = Epetra.Vector(Map);
    RHS = Epetra.Vector(Map);
    ExactSolution = Epetra.Vector(Map); ExactSolution.Random();
    Matrix.Apply(ExactSolution, RHS);
    LHS.PutScalar(0.0);

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

  step = sys.argv[2]

  if step == 'analyze':
    counter = 0;
    donothing = 0;
  elif step == 'iterative':
    counter = int(config['COUNTER']);
    List['solver'] = config['SOLVER'];
    List['iters'] = int(config['ITERS']);
    List['tol'] = float(config['TOL']);
  elif step == 'direct':
    counter = int(config['COUNTER']);
    List['solver'] = config['SOLVER'];

  study_parameter = {}

  for i in range(1,counter-1):
    type_and_name = config['__PyTrilinos__name' + str(i)];
    type = type_and_name[0:type_and_name.find(':')]
    name = type_and_name[type_and_name.find(':')+1:]
    value = config['__PyTrilinos__value' + str(i)];

    if len(value.split('; ')) > 1:
      study_parameter[type_and_name] = value
    else:  
      set_type(List, name, type, value)

  # Construct the problem =====================================================

  print '>>> ', config['PROBLEM_ID']
  problemIDs = config['PROBLEM_ID'].split(":");
  
  for problemID in problemIDs:
    if problemID == '':
      continue;
    print '### PROBLEM ', problemID, ' ###'
    print problemID
    (Map, Matrix, LHS, RHS, ExactSolution) = generator(problemID, comm);

    if step == 'analyze':
      analyze(Map, Matrix, LHS, RHS, ExactSolution);
    elif step == 'iterative':
      if len(study_parameter.keys()) == 0:
        iterative(Map, Matrix, LHS, RHS, ExactSolution, List);
      elif len(study_parameter.keys()) == 1:
        type_and_name = study_parameter.keys()[0]

        type = type_and_name[0:type_and_name.find(':')]
        name = type_and_name[type_and_name.find(':')+1:]
        for value in study_parameter[type_and_name].split('; '):
          print 'KEY = ', name, ' VALUE = ', value
          set_type(List, name, type, value)
          print 'ATTENTION: ZERO INITIAL SOL'
          LHS.PutScalar(0.0)
          iterative(Map, Matrix, LHS, RHS, ExactSolution, List);
      elif len(study_parameter.keys()) == 2:
        type_and_name0 = study_parameter.keys()[0]
        type_and_name1 = study_parameter.keys()[1]

        type0 = type_and_name0[0:type_and_name0.find(':')]
        name0 = type_and_name0[type_and_name0.find(':')+1:]

        type1 = type_and_name1[0:type_and_name1.find(':')]
        name1 = type_and_name1[type_and_name1.find(':')+1:]

        for value0 in study_parameter[type_and_name0].split('; '):
          for value1 in study_parameter[type_and_name1].split('; '):
            print 'KEY_0 = ', name0, ' VALUE_0 = ', value0
            print 'KEY_1 = ', name1, ' VALUE_1 = ', value1

            set_type(List, name0, type0, value0)
            set_type(List, name1, type1, value1)

            print 'ATTENTION: ZERO INITIAL SOL'
            LHS.PutScalar(0.0)
            iterative(Map, Matrix, LHS, RHS, ExactSolution, List);
      else:
        print 'THIS IS ONLY ALLOWED WITH ONE or TWO PARAMETER'
    elif step == 'direct':
      direct(Map, Matrix, LHS, RHS, ExactSolution, List);

# -------------------------------------------------------------------------
if __name__ == "__main__":
  main()
