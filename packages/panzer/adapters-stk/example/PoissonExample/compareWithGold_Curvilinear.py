import math
import sys
import re

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

  tol = 1e-9

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

  for i in range(len(q2gold_errors)):
    if abs(q2gold_errors[i] - q2_errors[i]) > tol:
      print( 'Poisson example using a Q2 mesh does not pass regression check:',abs(q2gold_errors[i] - q2_errors[i]),"must be less than",tol )
      raise 'Exception'

    if abs(q1gold_errors[i] - q1_errors[i]) > tol:
      print( 'Poisson example using a Q1 mesh does not pass regression check:',abs(q1gold_errors[i] - q1_errors[i]),"must be less than",tol )
      raise 'Exception'

    if q2_errors[i] > q1_errors[i]:
      print ( 'Errors are larger when using Q2 mesh!' )
      raise 'Exception'

  print( 'Poisson example using a Q2 mesh passes regression check' )
  print( 'Poisson example using a Q1 mesh passes regression check' )
  print ( 'Errors are smaller when using Q2 mesh, as expected')

except Exception as e:
  print("Test Failed: ")
  print(e)
else:
  print("Test Passed")
