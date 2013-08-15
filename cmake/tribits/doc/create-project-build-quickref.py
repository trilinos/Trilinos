#!/usr/bin/env python

usageHelp = r"""create-project-build-quickref.py [options]

Generate a TriBITS project-specific build quick reference file.

ToDo: Finish documentation!

"""

import sys
import os
import stat
import subprocess
import commands

from optparse import OptionParser

#
# A) Set up basic paths and import modules
#

thisFilePath = __file__
thisFileRealAbsBasePath = \
  os.path.dirname(os.path.abspath(os.path.realpath(thisFilePath)))
#print "thisFileRealAbsBasePath = '"+thisFileRealAbsBasePath+"'"

tribitsBaseDir = os.path.abspath(os.path.join(thisFileRealAbsBasePath, '..'))

sys.path.append(os.path.join(tribitsBaseDir, 'python'))
#print "sys.path =", sys.path

from GeneralScriptSupport import *
from CheckinTest import *
import GenerateDocUtilsOutput


#
# B) Read in the commandline options
#
  
clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--project-base-dir", dest="projectBaseDir", type="string",
  default="",
  help="Base directory of the TriBITS project.  If not set, is set by default" \
    " to <tribitsBaseDir>/../.. which assumes that Tribits is under" \
    " cmake/tribits."
  )
   
clp.add_option(
  "--min-cmake-version", dest="minCMakeVersion", type="string",
  help="Minimum version of CMake needed for given project" \
    " [Default '2.8.1']",
  default="2.8.1" )

GenerateDocUtilsOutput.addCmndLineOptions(clp)

(options, args) = clp.parse_args(sys.argv)

#
# C) Get information for everything
#

if options.projectBaseDir:
  projectBaseDir = options.projectBaseDir
else:
  projectBaseDir = \
    os.path.abspath(os.path.join(tribitsBaseDir, "../.."))
#print "projectBaseDir =", projectBaseDir

projectName = getProjectName(projectBaseDir)
#print "projectName =", projectName

if options.fileBase:
  options.fileBase = options.fileBase
else:
  options.fileBase = os.path.join(projectBaseDir,
    projectName+"BuildQuickRef")
#print "options.fileBase =", options.fileBase

#
# D) Read in standard body and make substitution
#

tribitsBuildQuickRefBodyFile = \
  os.path.join(thisFileRealAbsBasePath, "TribitsBuildQuickRefBody.rst")
#print "tribitsBuildQuickRefBodyFile =", tribitsBuildQuickRefBodyFile

tribitsBuildQuickRefBodyStr = \
  readStrFromFile(tribitsBuildQuickRefBodyFile)

substitutedTribitsBuildQuickRefBodyStr = tribitsBuildQuickRefBodyStr
substitutedTribitsBuildQuickRefBodyStr = \
    substitutedTribitsBuildQuickRefBodyStr.replace("<Project>",
      projectName)
substitutedTribitsBuildQuickRefBodyStr = \
    substitutedTribitsBuildQuickRefBodyStr.replace("<MinCMakeVer>",
      options.minCMakeVersion)

#
# E) Generate the output files
#

projectBuildQuickRefTemplateFile = \
  os.path.join(projectBaseDir, "cmake", \
    projectName+"BuildQuickRefTemplate.rst")
#print "projectBuildQuickRefTemplateFile =", projectBuildQuickRefTemplateFile

projectBuildQuickRefTemplateStr = \
  readStrFromFile(projectBuildQuickRefTemplateFile)

projectBuildQuickRefStr = \
  projectBuildQuickRefTemplateStr \
  + "\n\n" \
  + substitutedTribitsBuildQuickRefBodyStr

fileBaseName = os.path.basename(options.fileBase)

outputRstFile = options.fileBase+".rst"
print "Writing rst file ..."
GenerateDocUtilsOutput.openWriteFilePermissions(outputRstFile)
writeStrToFile(outputRstFile, projectBuildQuickRefStr)
GenerateDocUtilsOutput.setGeneratedFilePermissions(outputRstFile)

GenerateDocUtilsOutput.generateDocutilsOuputFiles(options)
