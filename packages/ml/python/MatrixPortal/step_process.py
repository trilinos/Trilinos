#! /usr/bin/env python
import os
import sys
import re
import string
import Numeric

from PyTrilinos import Epetra, EpetraExt, AztecOO, ML, Galeri, IFPACK

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
      p2 = string.split(parts[i], '+')
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
def perform_analysis(Label, Map, Matrix, LHS, RHS, ExactSolution):
  print "<p><p><div class=\"outputBox\"><pre>";
  print "<b><font color=red>Problem Label = ", Label, "</font></b>";
  print "<b><font color=red>Operation = matrix analysis </font></b>";
  IFPACK.AnalyzeMatrix(Matrix);
  IFPACK.AnalyzeMatrixElements(Matrix);
  print "&nbsp;<pre></div>";

# -------------------------------------------------------------------------
def perform_IFPACK(What, Label, Map, Matrix, LHS, RHS, ExactSolution, List):
  print "<p><p><div class=\"outputBox\"><pre>";
  print "<b><font color=red>Problem Label = ", Label, "</font></b>";
  print "<b><font color=red>Operation = ", What, "</font></b>";
  Factory = IFPACK.Factory()
  if What == "Jacobi":
    List['relaxation: type'] = "Jacobi";
    Prec = Factory.Create("point relaxation stand-alone", Matrix)
  elif What == "Gauss-Seidel":
    List['relaxation: type'] = "Gauss-Seidel";
    Prec = Factory.Create("point relaxation stand-alone", Matrix)
  elif What == "symmetric Gauss-Seidel":
    List['relaxation: type'] = "symmetric Gauss-Seidel";
    Prec = Factory.Create("point relaxation stand-alone", Matrix)
  elif What == "IC":
    Prec = Factory.Create("IC stand-alone", Matrix);
  elif What == "ICT":
    Prec = Factory.Create("ICT stand-alone", Matrix);
  elif What == "ILU":
    Prec = Factory.Create("ILU stand-alone", Matrix);
  elif What == "ILUT":
    Prec = Factory.Create("ILUT stand-alone", Matrix);

  Prec.SetParameters(List)
  Prec.Initialize()
  Prec.Compute()

  RHS.Random();
  LHS.PutScalar(0.0);
  
  Solver = AztecOO.AztecOO(Matrix, LHS, RHS)
  Solver.SetPrecOperator(Prec)
  if (List['az_solver'] == "AZ_gmres"):
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_gmres);
  elif List['az_solver'] == "AZ_cg":
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg);
  elif List['az_solver'] == "AZ_cg_condnum":
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg_condnum);
  elif List['az_solver'] == "AZ_gmres_condnum":
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_gmres_condnum);
  elif List['az_solver'] == "AZ_bicgstab":
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_bicgstab);
  elif List['az_solver'] == "AZ_tfqmr":
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_tfqmr);
  else:
    print "Solver type not correct, ", List['az_solver']
  Solver.SetAztecOption(AztecOO.AZ_output, 16);
  err = Solver.Iterate(List['iters'], List['tol']) 

  del Prec;
  print "&nbsp;<pre></div>";
  return(Solver.NumIters())
  
# -------------------------------------------------------------------------
def perform_ml(Label, Map, Matrix, LHS, RHS, ExactSolution, List):
  print "<p><p><div class=\"outputBox\"><pre>";
  print "<b><font color=red>Problem Label = ", Label, "</font></b>";
  print "<b><font color=red>Operation = multilevel preconditioner </font></b>";
  Prec = ML.MultiLevelPreconditioner(Matrix, False);
  Prec.SetParameterList(List);
  Prec.ComputePreconditioner();

  RHS.Random();
  LHS.PutScalar(0.0);
  
  Solver = AztecOO.AztecOO(Matrix, LHS, RHS)
  Solver.SetPrecOperator(Prec)
  if (List['az_solver'] == "AZ_gmres"):
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_gmres);
  elif List['az_solver'] == "AZ_cg":
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg);
  elif List['az_solver'] == "AZ_cg_condnum":
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg_condnum);
  elif List['az_solver'] == "AZ_gmres_condnum":
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_gmres_condnum);
  elif List['az_solver'] == "AZ_bicgstab":
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_bicgstab);
  elif List['az_solver'] == "AZ_tfqmr":
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_tfqmr);
  else:
    print "Solver type not correct, ", List['az_solver']
  Solver.SetAztecOption(AztecOO.AZ_output, 16);
  err = Solver.Iterate(List['iters'], List['tol']) 

  del Prec;
  print "&nbsp;<pre></div>";
  return(Solver.NumIters())
  
# -------------------------------------------------------------------------
def direct(Map, Matrix, LHS, RHS, ExactSolution, List):
  List['PrintStatus'] = True;
  List['PrintTiming'] = True;

  Problem = Epetra.LinearProblem(Matrix, LHS, RHS);
  Factory = Amesos.Factory();
  Solver = Factory.Create(List['direct_solver'], Problem);
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
  List = {}
  for l in file.readlines():
    d = string.split(string.strip(l), '=')
    if string.strip(d[0]) == "ProblemIDs":
      ProblemIDs = d[1];
      continue;
    what = d[0][0];
    name = string.strip(d[0][2:]);
    val = d[1];
    if what == "i":
      List[name] = int(val);
    elif what == "b":
      List[name] = bool(val);
    elif what == "d":
      List[name] = float(val);
    elif what == "s":
      List[name] = string.strip(val);

  # Construct the problem =====================================================

  ProblemIDs = ProblemIDs[1:].split(":");
  
  for FullProblemID in ProblemIDs:
    if FullProblemID == '':
      continue;

    pos = FullProblemID.find('@');
    Label = FullProblemID[0:pos - 1];
    problemID = FullProblemID[pos + 1:];
    (Map, Matrix, LHS, RHS, ExactSolution) = generator(problemID, comm);

    if List['perform_analysis']:
      perform_analysis(Label, Map, Matrix, LHS, RHS, ExactSolution);
    
    if List['perform_jacobi']:
      phi = perform_IFPACK("Jacobi", Label, Map, Matrix, LHS, RHS, ExactSolution, List);
    
    if List['perform_gs']:
      phi = perform_IFPACK("Gauss-Seidel", Label, Map, Matrix, LHS, RHS, ExactSolution, List);
    
    if List['perform_sgs']:
      phi = perform_IFPACK("symmetric Gauss-Seidel", Label, Map, Matrix, LHS, RHS, ExactSolution, List);
    
    if List['perform_ic']:
      phi = perform_IFPACK("IC", Label, Map, Matrix, LHS, RHS, ExactSolution, List);
    
    if False & List['perform_ict']: # FIXME
      phi = perform_IFPACK("ICT", Label, Map, Matrix, LHS, RHS, ExactSolution, List);
    
    if List['perform_ilu']:
      phi = perform_IFPACK("ILU", Label, Map, Matrix, LHS, RHS, ExactSolution, List);
    
    if List['perform_ilut']:
      phi = perform_IFPACK("ILUT", Label, Map, Matrix, LHS, RHS, ExactSolution, List);
    
    if List['perform_ml']:
      phi = perform_ml(Label, Map, Matrix, LHS, RHS, ExactSolution, List);

# -------------------------------------------------------------------------
if __name__ == "__main__":
  main()
