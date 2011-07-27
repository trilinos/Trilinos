#!/usr/bin/env python

import sys
import os
import re
import traceback


scriptsDir = os.path.abspath(os.path.dirname(sys.argv[0]))+"/../../cmake/python"
sys.path.insert(0, scriptsDir)

from GeneralScriptSupport import *


#
# Read in the command-line arguments
#

usageHelp = r"""update-copyright-header.py [OPTIONS]

This script updates the copyright header in a given file.  The logic in looks
for an existing copyright header in the file and replaces it with the provided
copyright header.  A copyright header is determined by looking for beginning
@HEADER and ending @HEADER lines.  If no such comment blocks are found, then
the given copyright header is inserted as the very first line in the file
enclosed in '/*' and '*/' lines so as to work also for C files.

To replace the copyright header for all of the source files for a complete
package do:

  $ cd packages/SOME_PACKAGE
  $ find . -name "*pp" -exec \
      ../../commonTools/refactoring/update-copyright-header.py  \
        --copyright-header=$PWD/Copyright.txt --file={} \;

For example, the file packages/SOME_PACKAGE/Copyright.txt should look like:

  $ cat packages/teuchos/Copyright.txt
    // @HEADER
    // ***********************************************************************
    // 
    //                    Teuchos: Common Tools Package
    //                 Copyright (2004) Sandia Corporation
    //
    // ...
    // 
    // ***********************************************************************
    // @HEADER

See packages/teuchos/Copyright.txt for an example.
"""

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--copyright-header", dest="copyrightHeader", type="string", default="",
  help="File containing the copyright header the be used." )

clp.add_option(
  "--file", dest="fileName", type="string", default="",
  help="File to have the copyright header replaced in." )

clp.add_option(
  "--script-mode", dest="scriptMode", action="store_true", default=False,
  help="Enable script mode. This prevents the addition of c style comment " +
       "markers at the begining and end of the copyright notice" )

(options, args) = clp.parse_args()


#
# Check/adjust the input arguments
#

if not options.copyrightHeader:
  raise Exception("Error, must set --copyrightHeader")

if not options.fileName:
  raise Exception("Error, must set --fileName")

if options.scriptMode == False:
  commentBegin = "/*" + os.linesep
  commentEnd   = "*/" + os.linesep
else:
  commentBegin = ""
  commentEnd   = ""


#
# Executble code
#

# A) Read in the standard copyright header

copyrightHeaderStr = open(options.copyrightHeader, 'r').read()
#print "copyrightHeaderStr:\n----------\n" + copyrightHeaderStr+"-----------\n"

# B) Read in the given file line by line looking to replace the copyright header

reHeaderLine = re.compile(r".*@HEADER.*")

fileLines = open(options.fileName, 'r').readlines()

newFileStr = ""

foundExistingHeaderBlock = False
inHeaderBlock = False

for line in fileLines:
  #print "line: '"+line+"'"
  if reHeaderLine.match(line) and not inHeaderBlock:
    #print "Found first @HEADER!"
    if foundExistingHeaderBlock:
      raise Exception("Error, the file '"+options.fileName+"'" + \
        " contains more than one header block!");
    foundExistingHeaderBlock = True
    inHeaderBlock = True
    newFileStr += copyrightHeaderStr
  elif reHeaderLine.match(line) and inHeaderBlock:
    #print "Found last @HEADER!"
    inHeaderBlock = False
  elif inHeaderBlock:
    #print "Skipping line in existing header block!"
    None
  else:
    #print "Augmenting line!"
    newFileStr += line

# C) If an existing header was never found, then add one at the top of the
# file.

if not foundExistingHeaderBlock:
  newFileStr = commentBegin + copyrightHeaderStr + commentEnd + os.linesep + \
               newFileStr


# D) Write the new file

#print "\n\nnewFileStr:\n----------\n"+newFileStr+"-----------\n"
open(options.fileName, 'w').write(newFileStr)
