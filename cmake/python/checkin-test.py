#!/usr/bin/env python

"""
Script for doing checkin testing of Trilinos.  Please run checkin-test.py -h
for details
"""


#
# Import commands
#


import sys
import os
import traceback

from GeneralScriptSupport import *

#
# Read in the commandline arguments
#

#print "sys.argv:", sys.argv

# Create a deep copy of the commandline arguments
cmndLineArgs = []
cmndLineArgs.extend(sys.argv)

# See if the help option is set or not
helpOpt = len( set(cmndLineArgs) & set(("--help", "-h")) ) > 0

# See if --show-defaults was set or not
showDefaultsOpt = len( set(cmndLineArgs) & set(("--show-defaults", "dummy")) ) > 0

#
# Forward the options but tee the output
#

scriptsDir = getScriptBaseDir()

if (not helpOpt) and (not showDefaultsOpt):
  logFileName = "checkin-test.out"
else:
  logFileName = ""

cmnd = scriptsDir+"/checkin-test-impl.py " + requoteCmndLineArgs(sys.argv[1:])

if logFileName:
  cmnd = cmnd + " 2>&1 | tee "+logFileName

rtn = echoRunSysCmnd(cmnd, throwExcept=False)

exit(rtn)
