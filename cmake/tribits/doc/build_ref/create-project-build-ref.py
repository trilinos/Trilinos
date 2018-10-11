#!/usr/bin/env python

usageHelp = r"""create-project-build-ref.py [options]

Generate a TriBITS project-specific build quick reference guide, typically
named <projectName>BuildReference.rst and output formats *.html and *.pdf (by
default).  For the average TriBITS project-specific build quick reference
file, you just run this as::

  $ cd <projectBaseDir>
  $ ./cmake/tribits/doc/create-project-build-ref.py

(or run create-project-build-ref.py from where TriBITS is located for the
project).  This will read the project's name from the file
<projectBaseDir>/ProjectName.cmake and will take the template file:

   <projectBaseDir>/cmake/<projectName>BuildReferenceTemplate.rst

and the general TribitsBuildReferenceBody.rst input file and generate the
files:

   <projectBaseDir>/
     <projectName>BuildReference.rst
     <projectName>BuildReference.html
     <projectName>BuildReference.pdf

However, this script can be run to take any input template *.rst file
<projectTemplateFile> and any arbitrary project name and use it to generate
the specified output files.  For an example of this, see the script:

  tribits/doc/build_quick_ref/create-build-ref.sh

that creates a the files TribitsBuildReference.[rst,html,pdf]

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

tribitsBaseDir = os.path.abspath(os.path.join(thisFileRealAbsBasePath, '../..'))
pythonUtilsDir = os.path.join(tribitsBaseDir, 'python_utils')
ciSupportDir = os.path.join(tribitsBaseDir, 'ci_support')

sys.path = [ciSupportDir, pythonUtilsDir] + sys.path
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
  help="Base directory of the TriBITS project <projectBaseDir>.  If not set, is set"+\
    " by default <tribitsBaseDir>/../.. which assumes that Tribits is under" \
    " cmake/tribits in the current project.  This option is not used if --project-name,"+\
    " --project-template-file, and --file-base are all overridden."
  )

clp.add_option(
  "--project-name", dest="projectName", type="string",
  default="",
  help="If specified, then this will be ued for the project name (<projectName>."+\
    "  Otherwise, the project name is read from <projectBaseDir>/ProjectName.cmake instead."
  )

clp.add_option(
  "--project-template-file", dest="projectTemplateFile", type="string",
  default="",
  help="The project-specific template file used to generate the overall *.rst file."+\
    "   If not specified, then <projectBaseDir>/cmake/<projectName>BuildReferenceTemplate.rst"+\
    " will be used instead."
  )
   
clp.add_option(
  "--min-cmake-version", dest="minCMakeVersion", type="string",
  help="Minimum version of CMake needed for given project" \
    " [Default '3.10.0']",
  default="3.10.0" )

GenerateDocUtilsOutput.addCmndLineOptions(clp)

(options, args) = clp.parse_args(sys.argv)

#
# C) Get information for everything
#

if not options.projectBaseDir:
  options.projectBaseDir = \
    os.path.abspath(os.path.join(tribitsBaseDir, "../.."))
#print "projectBaseDir =", projectBaseDir

if not options.projectName:
  options.projectName = getProjectName(options.projectBaseDir)
#print "projectName =", projectName

if not options.projectTemplateFile:
  options.projectTemplateFile = \
    os.path.join(options.projectBaseDir, "cmake", \
      options.projectName+"BuildReferenceTemplate.rst")
#print "projectBuildReferenceTemplateFile =", projectBuildReferenceTemplateFile

if not options.fileBase:
  options.fileBase = os.path.join(options.projectBaseDir,
    options.projectName+"BuildReference")
#print "options.fileBase =", options.fileBase

#
# D) Read in standard body and make substitution
#

tribitsBuildReferenceBodyFile = \
  os.path.join(thisFileRealAbsBasePath, "TribitsBuildReferenceBody.rst")
#print "tribitsBuildReferenceBodyFile =", tribitsBuildReferenceBodyFile

tribitsBuildReferenceBodyStr = \
  readStrFromFile(tribitsBuildReferenceBodyFile)

substitutedTribitsBuildReferenceBodyStr = tribitsBuildReferenceBodyStr
substitutedTribitsBuildReferenceBodyStr = \
    substitutedTribitsBuildReferenceBodyStr.replace("<Project>",
      options.projectName)
substitutedTribitsBuildReferenceBodyStr = \
    substitutedTribitsBuildReferenceBodyStr.replace("<MinCMakeVer>",
      options.minCMakeVersion)

#
# E) Generate the output files
#

projectBuildReferenceTemplateStr = \
  readStrFromFile(options.projectTemplateFile)

projectBuildReferenceStr = \
  projectBuildReferenceTemplateStr \
  + "\n\n" \
  + substitutedTribitsBuildReferenceBodyStr

outputRstFile = options.fileBase+".rst"
print "Writing rst file ..."
GenerateDocUtilsOutput.openWriteFilePermissions(outputRstFile)
writeStrToFile(outputRstFile, projectBuildReferenceStr)
GenerateDocUtilsOutput.setGeneratedFilePermissions(outputRstFile)

GenerateDocUtilsOutput.generateDocutilsOuputFiles(options)
