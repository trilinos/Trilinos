#! /usr/bin/env python
import sys
import re
import Numeric

try:
  import setpath
  import Epetra
  import EpetraExt
  import AztecOO
  import ML
except ImportError:
  from PyTrilinos import Epetra, EpetraExt, AztecOO, ML
  print "Using installed versions of Epetra, Galeri, AztecOO, ML"

###############################################################################
# main()
###############################################################################

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
    #print key+": "+config[key]
  
  # Construct the problem =====================================================

  #=== Comm ===#
  comm = Epetra.PyComm()
  
  #=== Map ===#
  
  # These two lines work correctly...
  n = 19881
  map = Epetra.Map(n, 0, comm)
  
  # ...this doesn't...but we'd like to create the map directly from the file...
  #(ierr, map) = EpetraExt.MatrixMarketFileToMap(config['A_FILE'], comm)
  
  #=== Matrix ===#
  (ierr, matrix) = EpetraExt.MatrixMarketFileToCrsMatrix(config['A_FILE'], map)
  
  #=== LHS ===#
  if config.has_key('X_FILE') and config['X_FILE']:
    file = open(config['X_FILE'], 'r')
    fileContents = file.read()  
    pattern = re.compile(r"\s*\d+\s+\d+\s+((?:-|\d|.)+)\s*$", re.M)
    values = pattern.findall(fileContents)
    values = [float(v) for v in values]
    lhs = Epetra.Vector(map, Numeric.array(values))
  else:
    lhs = Epetra.Vector(map); lhs.Random()
  
  #=== RHS ===#
  if config.has_key('B_FILE') and config['B_FILE']:
    file = open(config['B_FILE'], 'r')
    fileContents = file.read()  
    pattern = re.compile(r"\s*\d+\s+\d+\s+((?:-|\d|.)+)\s*$", re.M)
    values = pattern.findall(fileContents)
    values = [float(v) for v in values]
    rhs = Epetra.Vector(map, Numeric.array(values))
  else:
    rhs = Epetra.Vector(map); rhs.PutScalar(0.0)
  
  # Solve the problem =========================================================

  # sets up the parameters for ML using a python dictionary
  mlList = {}
  
  mlList['output'] = 10
  
  if config.has_key('SMOOTHER_TYPE') and config['SMOOTHER_TYPE']:
    mlList['smoother: type'] = config['SMOOTHER_TYPE']
    
  if config.has_key('SMOOTHER_SWEEPS') and config['SMOOTHER_SWEEPS']:
    mlList['smoother: sweeps'] = int(config['SMOOTHER_SWEEPS'])
    
  if config.has_key('SMOOTHER_MLS_POLYNOMIAL_ORDER') and config['SMOOTHER_MLS_POLYNOMIAL_ORDER']:
    mlList['smoother: MLS plynomial order'] = int(config['SMOOTHER_MLS_POLYNOMIAL_ORDER'])
    
  if config.has_key('COARSE_TYPE') and config['COARSE_TYPE']:
    mlList['coarse: type'] = config['COARSE_TYPE']
    
  if config.has_key('COARSE_SWEEPS') and config['COARSE_SWEEPS']:
    mlList['coarse: sweeps'] = int(config['COARSE_SWEEPS'])
    
  if config.has_key('COARSE_MLS_POLYNOMIAL_ORDER') and config['COARSE_MLS_POLYNOMIAL_ORDER']:
    mlList['coarse: MLS plynomial order'] = int(config['COARSE_MLS_POLYNOMIAL_ORDER'])
    
  if config.has_key('MAX_LEVELS') and config['MAX_LEVELS']:
    mlList['max levels'] = int(config['MAX_LEVELS'])
    
  if config.has_key('AGGREGATION_TYPE') and config['AGGREGATION_TYPE']:
    mlList['aggregation: type'] = config['AGGREGATION_TYPE']
    
  if config.has_key('AGGREGATION_DAMPING_FACTOR') and config['AGGREGATION_DAMPING_FACTOR']:
    mlList['aggregation: damping factor'] = float(config['AGGREGATION_DAMPING_FACTOR'])

  # creates the preconditioner and computes it
  prec = ML.MultiLevelPreconditioner(matrix, False)
  prec.SetParameterList(mlList)
  prec.ComputePreconditioner()

  # sets up the solver, specifies Prec as preconditioner, and solves using CG.
  solver = AztecOO.AztecOO(matrix, lhs, rhs)
  solver.SetPrecOperator(prec)
  solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg);
  solver.SetAztecOption(AztecOO.AZ_output, 16);
  solver.Iterate(1550, 1e-5)

###############################################################################

if __name__ == "__main__":
  main()
