import math
import sys
import re
import numpy as np

try:
  if not len(sys.argv)==5:
    print('Expected four files, two gold and the two test outputs')
    raise Exception("Failure: Command line - <cmd> [test_Q2_mesh.gold] [test_Q1_mesh.gold] [test_Q2_mesh.out] [test_Q1_mesh.out]")

  q2gold = sys.argv[1]
  q1gold = sys.argv[2]
  q2test = sys.argv[3]
  q1test = sys.argv[4]

  ErrorPat = re.compile('.*Error = (.*)$')

  q2gold_errors = []
  q1gold_errors = []
  q2_errors = []
  q1_errors = []

  with open(q2gold) as f:
    for line in f:
      match = ErrorPat.match(line)
      if match!=None:
        q2gold_errors.append(float(match.group(1)))
  with open(q1gold) as f:
    for line in f:
      match = ErrorPat.match(line)
      if match!=None:
        q1gold_errors.append(float(match.group(1)))
  with open(q2test) as f:
    for line in f:
      match = ErrorPat.match(line)
      if match!=None:
        q2_errors.append(float(match.group(1)))
  with open(q1test) as f:
    for line in f:
      match = ErrorPat.match(line)
      if match!=None:
        q1_errors.append(float(match.group(1)))

  q2gold_errors = np.asarray(q2gold_errors)
  q1gold_errors = np.asarray(q1gold_errors)
  q2_errors = np.asarray(q2_errors)
  q1_errors = np.asarray(q1_errors)

  if (np.any( q2gold_errors != q2_errors ) ):
    print( 'Poisson example using a Q2 mesh does not pass regression check' )
    raise 'Exception'
  else:
    print( 'Poisson example using a Q2 mesh passes regression check' )

  if (np.any( q1gold_errors != q1_errors ) ):
    print( 'Poisson example using a Q1 mesh does not pass regression check' )
    raise 'Exception'
  else:
    print( 'Poisson example using a Q1 mesh passes regression check' )

  if (np.any( q2_errors > q1_errors ) ):
    print ( 'Errors are larger when using Q2 mesh!' )
    raise 'Exception'
  else:
    print ( 'Errors are smaller when using Q2 mesh, as expected')

except Exception as e:
  print("Test Failed: ")
  print(e)
else:
  print("Test Passed")
