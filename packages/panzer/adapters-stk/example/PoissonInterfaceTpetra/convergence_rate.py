import sys
import re

try:
  if not len(sys.argv)>3:
    print('Expected at least three solution files with the line \"Error = \%f\" in them')
    raise Exception("Failure: Command line - <cmd> [filename-prefix] [coarsest mesh res] [second mesh res] [third mesh res] ...")

  file_prefix = sys.argv[1]

  mesh_res = []
  for i in range(2,len(sys.argv)):
    mesh_res += [int(sys.argv[i])]

  max_res =  max(mesh_res)
  
  errorPat = re.compile('.*Error = (.*)$')

  error_values = []
  for res in mesh_res:
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
      errorMatch = errorPat.match(line)
      if errorMatch==None:
        continue

      error = float(errorMatch.group(1))
      error_values += [error]
  # end for res

  rate = error_values[-1]/error_values[-2]

  # the convergence rate should be 0.25 (we are halfing the mesh size
  # each refinement step). Though its not exact due to using tets and
  # other reasons that lead to imperfect convergence rates. But its close
  # enough to give a quadratic convergence order.
  if (rate-0.25)/0.25 > 0.05:
    raise Exception('Convergence rate of %f, is not within 5 percent bounds of 0.25: FAILED' % rate)

  print('Convergence rate of %f is within 5 percent of 0.25: PASSED' % rate)

except Exception as e:
  print("Test Failed: ")
  print(e)
else:
  print("Test Passed")
