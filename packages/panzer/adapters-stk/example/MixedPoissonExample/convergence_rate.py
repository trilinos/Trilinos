import math
import sys
import re

try:
  if not len(sys.argv)>4:
    print('Expected at least three solution files with the line \"Error = \%f\" in them')
    raise Exception("Failure: Command line - <cmd> [basis-order] [filename-prefix] [coarsest mesh res] [second mesh res] [third mesh res] ...")

  basis_order = float(sys.argv[1])
  file_prefix = sys.argv[2]

  mesh_res = []
  for i in range(3,len(sys.argv)):
    mesh_res += [int(sys.argv[i])]

  max_res =  max(mesh_res)

  #orderPat   = re.compile('.*Basis Order = (.*)$')  
  ErrorPat = re.compile('.*Error = (.*)$')

  error_value = 0
  prev_error_value = 0

  error_values = []
  rate_values = []

  for cnt, res in enumerate(mesh_res):
    if max_res>=100:
      filename = '%s%03d' % (file_prefix,res)
    elif max_res>=10:
      filename = '%s%02d' % (file_prefix,res)
    else:
      filename = '%s%d' % (file_prefix,res)

    print("opening \"" + filename +"\"")
    f = open(filename);

    # look for the "Error = ..." line
    for line in f:
      #orderMatch = orderPat.match(line)
      #if orderMatch!=None:
      #  basis_order = float(orderMatch.group(1))

      errorMatch = ErrorPat.match(line)
      if errorMatch!=None:
        error_value = float(errorMatch.group(1))
        error_values += [error_value]

    if cnt>0:
      prev_error = prev_error_value
      error = error_value
      rate = math.log(prev_error/error)/math.log(2)
      rate_values += [rate]
      diff = rate - basis_order
      # this is L2 error and converge with p+1 
      # 1st order element is too coarse for testing
      if not (basis_order == 1 and cnt <2):
        if rate > 0.90*(basis_order+1):
          print('%s: Convergence rate of %f is within 10 percent of %f: PASSED' % (filename, rate, basis_order+1)) 
        else:
          print('%s: Convergence rate of %f is not within 10 percent of %f: FAILED' % (filename, rate, basis_order+1)) 
          raise 'Exception' 

    prev_error_value = error_value

  # end for res

  print('Errors   : ', error_values)
  print('Conv rate   : ', rate_values)

except Exception as e:
  print("Test Failed: ")
  print(e)
else:
  print("Test Passed")
