#!/usr/bin/env python3

usageHelp = r"""generate-docutils-output.py --file-base=<FILEBASE> [options]

Create output files for DocUtils reStructuredText.

This simple script can be used to generate HTML, Latex, and PDF output for a
given <outputFileBase>.rst file.  What output is generated is controlled by
the options --generate-html, --generate-latex, and --generate-pdf.  This
requires that the Docutils commands be installed on the system as well as
latex (if generating PDF output).

NOTE: The right default programs for --generate-html, --generate-latex, and
--generate-pdf will try to be chosen for given system automatically.
"""

import sys
import os
from optparse import OptionParser

#
# A) Set up basic paths and import modules
#

thisFilePath = __file__
thisFileRealAbsBasePath = \
  os.path.dirname(os.path.abspath(os.path.realpath(thisFilePath)))
#print("thisFileRealAbsBasePath = '" + thisFileRealAbsBasePath + "'")

tribitsBaseDir = os.path.abspath(os.path.join(thisFileRealAbsBasePath, '..'))

sys.path.append(os.path.join(tribitsBaseDir, 'python'))
#print("sys.path = " + str(sys.path))


from GenerateDocUtilsOutput import *

#
# B) Read in the commandline options
#

clp = OptionParser(usage=usageHelp)

addCmndLineOptions(clp)

(options, args) = clp.parse_args(sys.argv)

if not options.fileBase:
  print("Error, --file-base=<fileBase> must be specified!")
  sys.exit(1)
  
#
# C) Generate the output files
#

generateDocutilsOuputFiles(options)
