#! /usr/bin/env python

import os
import sys

# set paths
import setpath

# Get vebosity from the command line
if (len(sys.argv) > 2):
    print "usage:", sys.argv[0], "[-v]"
    sys.exit(1)
if (len(sys.argv) == 2 and sys.argv[1] == "-v"):
    verbose = True
else:
   verbose = False
   
# run test
import testLOCA_Chan
result = testLOCA_Chan.runTests(verbose)

if (result == 0):
    print "Test was successful!"
else:
    print "Test Failed!"


sys.exit(result)
