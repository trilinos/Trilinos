#! /usr/bin/env python

print "Starting lapack/LOCA_python/testLOCA.exe"

import os
import sys

# set paths
import setpath

# run test
import testLOCA_Chan
result = testLOCA_Chan.main()

if (result == 0):
    print "Test was successful!"
else:
    print "Test Failed!"


sys.exit(result)
