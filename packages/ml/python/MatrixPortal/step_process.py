#! /usr/bin/env python
import os
import sys
import re
import string
import math

# -------------------------------------------------------------------------
# \author Marzio Sala, COLAB/ETHZ
#
# \date Last updated on Apr-06
# 
# To test this file on your local machine without any PHP/HTML, you
# can run the following command:
#
# $ python step_process.py step_text.txt
#
# which should run all the solvers with two test cases.
# -------------------------------------------------------------------------

from PyTrilinos import Epetra, EpetraExt, AztecOO, ML, Galeri, IFPACK

comm = Epetra.PyComm()

NumProcs = comm.NumProc()
MyPID = comm.MyPID()
count = 0
analysis_count = 0
  
# -------------------------------------------------------------------------
# checks on parameters and location of the HB and H5 files
# -------------------------------------------------------------------------
MAX_MATRIX_ROWS = 256 * 256
MAX_MATRIX_NONZEROS = MAX_MATRIX_ROWS * 5
MAX_ITERATIONS = 1550
MAX_KSPACE = 200
HB_REPOSITORY = "/home/chinella/Web/MatrixPortal/HBMatrices/"
H5_REPOSITORY = "/home/chinella/Web/MatrixPortal/H5Matrices/"

# -------------------------------------------------------------------------
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
def add_result(List, Label, ierr, iters, PrecTime, SolverTime, ConditionNumber):
  if List['evaluation: type'] == 'iterations':
    phi = iters;
  elif List['evaluation: type'] == 'condnum':
    phi = ConditionNumber
  elif List['evaluation: type'] == 'setup+sol':
    phi = PrecTime + SolverTime
  elif List['evaluation: type'] == 'setup':
    phi = PrecTime
  elif List['evaluation: type'] == 'sol':
    phi = SolverTime

  global count
  print "<p><font color=midnightblue>Evaluation phi = ",
  print "<input type=text size=10 name=phi_value_%d value=%f>" %(count, phi),
  print "Add to results with label: <input type=text size=20 name=phi_label_%d value=\"%s\"></font>" % (count, Label)
  count = count + 1

# -------------------------------------------------------------------------
def iterative_solver(List, Matrix, InputLHS, RHS, Prec):

  LHS = Epetra.MultiVector(InputLHS)
  
  Time = Epetra.Time(Matrix.Comm())

  hasConditionNumber = False;

  Solver = AztecOO.AztecOO(Matrix, LHS, RHS)
  Solver.SetPrecOperator(Prec)
  if (List['az_solver'] == "AZ_gmres"):
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_gmres);
  elif List['az_solver'] == "AZ_cg":
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg);
  elif List['az_solver'] == "AZ_cg_condnum":
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg_condnum);
    hasConditionNumber = True
  elif List['az_solver'] == "AZ_gmres_condnum":
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_gmres_condnum);
    hasConditionNumber = True
  elif List['az_solver'] == "AZ_bicgstab":
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_bicgstab);
  elif List['az_solver'] == "AZ_tfqmr":
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_tfqmr);
  elif List['az_solver'] == "AZ_cgs":
    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cgs);
  else:
    print "Solver type not correct, ", List['az_solver']

  Solver.SetAztecOption(AztecOO.AZ_output, 16);
  if List['iters'] < 0 | List['iters'] > MAX_ITERATIONS:
    print "Maximum number of iterations either negative of > ", MAX_ITERATIONS;
    raise("PARAMETER_ERROR");
  if List['tol'] < 1e-12:
    print "Tolerance is too small"
    raise("PARAMETER_ERROR");

  if List['az_kspace'] > 0 & List['az_kspace'] <= MAX_KSPACE:
    Solver.SetAztecOption(AztecOO.AZ_kspace, List['az_kspace']);
  else:
    print "Krylov space dimension either negative of >", MAX_KSPACE
    print "You have", List['az_kspace']
    raise("PARAMETER_ERROR");

  if List['az_output'] == "16":
    Solver.SetAztecOption(AztecOO.AZ_output, 16)
  elif List['az_output'] == "32":
    Solver.SetAztecOption(AztecOO.AZ_output, 32)
  elif List['az_output'] == "AZ_last":
    Solver.SetAztecOption(AztecOO.AZ_output, AztecOO.AZ_last)
  elif List['az_output'] == "AZ_none":
    Solver.SetAztecOption(AztecOO.AZ_output, AztecOO.AZ_none)
    
  err = Solver.Iterate(List['iters'], List['tol']) 

  if hasConditionNumber:
    ConditionNumber = Solver.GetStatus(AztecOO.AZ_condnum)
  else:
    ConditionNumber = 0.0;
  return (err, Solver.NumIters(), Time.ElapsedTime(), ConditionNumber)

# -------------------------------------------------------------------------
def generator(problemID, comm, List):
  GaleriList = {}
  if problemID[0:3] == "HB_":
    print "<p><p><div class=\"outputBox\"><pre>";
    FileName = HB_REPOSITORY + problemID[3:];
    Map, Matrix, LHS, RHS, ExactSolution = Galeri.ReadHB(FileName, comm);
    NullSpace = "not-set"
    print "</div>"

  elif problemID[0:3] == "XML_":
    print "TO BE DONE";

  else:
    parts = string.split(problemID, '_');
    ProblemType = parts[0];
    for i in range(1, len(parts)):
      p2 = string.split(parts[i], '!')
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
      mx = math.sqrt(NumProcs)
      if mx * mx == NumProcs:
        GaleriList['mx'] = int(mx)
        GaleriList['my'] = int(mx)
      else:
        GaleriList['mx'] = int(NumProcs)
        GaleriList['my'] = 1

    elif string.find(ProblemType, '3D') != -1:
      MapType = "Cartesian3D"
      mx = math.pow(NumProcs, 0.33334)
      if mx * mx * mx == NumProcs:
        GaleriList['mx'] = int(mx)
        GaleriList['my'] = int(mx)
        GaleriList['mz'] = int(mx)
      else:
        GaleriList['mx'] = int(NumProcs)
        GaleriList['my'] = 1
        GaleriList['mz'] = 1

    Map = Galeri.CreateMap(MapType, comm, GaleriList);

    Matrix = Galeri.CreateCrsMatrix(ProblemType, Map, GaleriList);

    LHS = Epetra.Vector(Map)
    RHS = Epetra.Vector(Map)
    ExactSolution = Epetra.Vector(Map)
    
    NullSpace = "not-set"

  # checks that the matrix is not too big

  if Map.NumGlobalElements() > MAX_MATRIX_ROWS:
    print "<b><font color=red>Sorry, the maximum matrix size is 20.000</font></b>"
    raise("PARAMETER_ERROR");

  if Matrix.NumGlobalNonzeros() > MAX_MATRIX_NONZEROS:
    print "<b><font color=red>Sorry, the maximum number of nonzeros is 250.000</font></b>"
    raise("PARAMETER_ERROR");

  # FIXME???
  if (List.has_key('solution') == True) & (List.has_key('starting_solution') == True) & (List.has_key('rhs') == True):

    if List['solution'] == "zero":
      ExactSolution.PutScalar(0.0)
    elif List['solution'] == "random":
      ExactSolution.Random()
    elif List['solution'] == "constant":
      ExactSolution.PutScalar(1.0)
    else:
      raise("PARAMETER_ERROR");

    if List['starting_solution'] == "zero":
      LHS.PutScalar(0.0)
    elif List['starting_solution'] == "random":
      LHS.Random()
    elif List['starting_solution'] == "constant":
      LHS.PutScalar(1.0)
    else:
      raise("PARAMETER_ERROR");
  
    if List['rhs'] == "zero":
      RHS.PutScalar(0.0)
    elif List['rhs'] == "random":
      RHS.Random()
    elif List['rhs'] == "constant":
      RHS.PutScalar(1.0)
    elif List['rhs'] == "matvec":
      Matrix.Apply(ExactSolution, RHS);
    else:
      raise("PARAMETER_ERROR");

  else:
    ExactSolution.Random()
    Matrix.Apply(ExactSolution, RHS)
    LHS.PutScalar(0.0)

  return(Map, Matrix, LHS, RHS, ExactSolution, NullSpace);

# -------------------------------------------------------------------------
def perform_analysis(Label, Map, Matrix, LHS, RHS, ExactSolution, List):
  if MyPID == 0:
    print "<p><p><div class=\"outputBox\"><pre>";
    print "<b><font color=red>Problem Label =", Label, "</font></b>";
    print "<b><font color=red>Operation = matrix analysis </font></b>";
  ImageBase = List['image_base'];
  TimeStamp = List['timestamp'];
  IFPACK.AnalyzeMatrix(Matrix, True);
  #IFPACK.AnalyzeMatrixElements(Matrix);
  global analysis_count
  PSImageBase = ImageBase + TimeStamp + str(analysis_count) + ".ps"
  PNGImageBase = ImageBase + TimeStamp + str(analysis_count) + ".png"

  IFPACK.PrintSparsity(Matrix, PSImageBase)
  import commands
  u=commands.getoutput('convert ' + PSImageBase + ' ' + PNGImageBase)
  if MyPID == 0:
    print '<center><p><img src=tmp/%s.png></center>' % (TimeStamp + str(analysis_count))
  analysis_count = analysis_count + 1
  
  if MyPID == 0:
    print "&nbsp;<pre></div>";

# -------------------------------------------------------------------------
def perform_IFPACK(What, Label, Map, Matrix, LHS, RHS, ExactSolution, List):
  if MyPID == 0:
    print "<p><p><div class=\"outputBox\"><pre>";
    print "<b><font color=midnightblue>Problem Label =", Label, "</font></b>";
    print "<b><font color=midnightblue>Operation = ", What, "</b>";
    print "Using", NumProcs, "processors."
    print "The AztecOO/IFPACK output follows."
    print "</font>";

  Time = Epetra.Time(Matrix.Comm())

  Factory = IFPACK.Factory()
  if What == "Chebyshev":
    Prec = Factory.Create("Chebyshev", Matrix)
  elif What == "Jacobi":
    List['relaxation: type'] = "Jacobi";
    Prec = Factory.Create("point relaxation stand-alone", Matrix)
  elif What == "Gauss-Seidel":
    List['relaxation: type'] = "Gauss-Seidel";
    Prec = Factory.Create("point relaxation stand-alone", Matrix)
  elif What == "symmetric Gauss-Seidel":
    List['relaxation: type'] = "symmetric Gauss-Seidel";
    Prec = Factory.Create("point relaxation stand-alone", Matrix)
  elif What == "IC":
    Prec = Factory.Create("IC", Matrix);
  elif What == "ICT":
    Prec = Factory.Create("ICT", Matrix);
  elif What == "ILU":
    Prec = Factory.Create("ILU", Matrix);
  elif What == "ILUT":
    Prec = Factory.Create("ILUT", Matrix);

  Prec.SetParameters(List)
  Prec.Initialize()
  Prec.Compute()

  PrecTime = Time.ElapsedTime()

  (ierr, iters, SolveTime, ConditionNumber) = iterative_solver(List, Matrix, LHS, RHS, Prec)
  del Prec;

  if MyPID == 0:
    add_result(List, Label + " " + What, ierr, iters, PrecTime, SolveTime, ConditionNumber)
    print "&nbsp;<pre></div>";
  
# -------------------------------------------------------------------------
def perform_ML(Label, Map, Matrix, LHS, RHS, ExactSolution, NullSpace, List):
  if MyPID == 0:
    print "<p><p><div class=\"outputBox\"><pre>";
    print "<b><font color=midnightblue>Problem Label =", Label, "</font></b>";
    print "<b><font color=midnightblue>Operation = multilevel preconditioner</b>";
    print "Using", NumProcs, "processors."
    print "The AztecOO/ML output follows."
    print "</font>";

  Time = Epetra.Time(Matrix.Comm())

  if NumProcs > 1:
    List['coarse: type'] = 'symmetric Gauss-Seidel';

  Prec = ML.MultiLevelPreconditioner(Matrix, False);
  if NullSpace == "not-set":
    Prec.SetParameterList(List)
  else:
    Prec.SetParameterListAndNullSpace(List, NullSpace)

  Prec.ComputePreconditioner()

  PrecTime = Time.ElapsedTime()

  (ierr, iters, SolveTime, ConditionNumber) = iterative_solver(List, Matrix, LHS, RHS, Prec)
  del Prec;

  if MyPID == 0:
    add_result(List, Label + " ML", ierr, iters, PrecTime, SolveTime, ConditionNumber)
    print "&nbsp;<pre></div>";
  
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
    d = string.split(string.strip(l), ':=')
    if d[0] == "":
      continue;
    if string.strip(d[0]) == "ProblemIDs":
      ProblemIDs = d[1];
      continue;
    what = d[0][0];
    name = string.strip(d[0][2:]);
    val = d[1];
    if (val != ""):
      if what == "i":
        List[name] = int(val)
      elif what == "b":
        List[name] = string.strip(val)
      elif what == "d":
        List[name] = float(val)
      elif what == "s":
        List[name] = string.strip(val)

  # Construct the problem =====================================================

  ProblemIDs = ProblemIDs[1:].split(":");
  
  for FullProblemID in ProblemIDs:
    if FullProblemID == '':
      continue;

    pos = FullProblemID.find('@');
    Label = FullProblemID[0:pos];
    problemID = FullProblemID[pos + 1:];
    (Map, Matrix, LHS, RHS, ExactSolution, NullSpace) = generator(problemID, comm, List);

    if List.has_key('perform_analysis'):
      if List['perform_analysis'] == "True":
        perform_analysis(Label, Map, Matrix, LHS, RHS, ExactSolution, List)
    
    if List.has_key('perform_cheby'):
      if  List['perform_cheby'] == "True":
        perform_IFPACK("Chebyshev", Label, Map, Matrix, LHS, RHS, ExactSolution, List);
    
    if List.has_key('perform_jacobi'):
      if  List['perform_jacobi'] == "True":
        perform_IFPACK("Jacobi", Label, Map, Matrix, LHS, RHS, ExactSolution, List);
    
    if List.has_key('perform_gs'):
      if  List['perform_gs'] == "True":
        perform_IFPACK("Gauss-Seidel", Label, Map, Matrix, LHS, RHS, ExactSolution, List);
    
    if List.has_key('perform_sgs'):
      if  List['perform_sgs'] == "True":
        perform_IFPACK("symmetric Gauss-Seidel", Label, Map, Matrix, LHS, RHS, ExactSolution, List);
    
    if List.has_key('perform_ic'):
      if  List['perform_ic'] == "True":
        perform_IFPACK("IC", Label, Map, Matrix, LHS, RHS, ExactSolution, List);
    
    if List.has_key('perform_ict'):
      if  List['perform_ict'] == "True": # FIXME
        perform_IFPACK("ICT", Label, Map, Matrix, LHS, RHS, ExactSolution, List);
    
    if List.has_key('perform_ilu') :
      if List['perform_ilu'] == "True":
        perform_IFPACK("ILU", Label, Map, Matrix, LHS, RHS, ExactSolution, List);
    
    if List.has_key('perform_ilut'):
      if List['perform_ilut'] == "True":
        perform_IFPACK("ILUT", Label, Map, Matrix, LHS, RHS, ExactSolution, List);
    
    if List.has_key('perform_ml'):
      if List['perform_ml'] == "True":
        perform_ML(Label, Map, Matrix, LHS, RHS, ExactSolution, NullSpace, List);

  if MyPID == 0:
    global count
    print '<input type=hidden name=phi_count value=%d>' % count

# -------------------------------------------------------------------------
if __name__ == "__main__":
  main()
