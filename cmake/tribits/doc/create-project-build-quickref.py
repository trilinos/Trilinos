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

#
# B) Define some helper functions
#

def openWriteFilePermissions(filePath):
  if os.path.exists(filePath):
    os.chmod(filePath, stat.S_IREAD | stat.S_IWRITE \
      | stat.S_IRGRP | stat.S_IWGRP)

def setGeneratedFilePermissions(filePath):
  os.chmod(filePath, stat.S_IREAD | stat.S_IRGRP | stat.S_IROTH)

def generateFile(filePath, generateCmnd, outFile=None, workingDir="", runTwice=False):
  openWriteFilePermissions(filePath)
  runSysCmnd(generateCmnd, outFile=outFile, workingDir=workingDir)
  if runTwice:
    runSysCmnd(generateCmnd, outFile=outFile, workingDir=workingDir)
  setGeneratedFilePermissions(filePath)


#
# C) Read in the commandline options
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
  "--output-file-base", dest="outputFileBase", type="string",
  default="",
  help="Project-specific build quickref file base name." \
    "  If not set, then set to <projectBaseDir>/<Project>BuildQuickRef"
  )
  
clp.add_option(
  "--generate-html", dest="generateHtml", type="string",
  help="Generate the HTML output file using provided script (i.e. rst2html)",
  default="rst2html.py" )
  
clp.add_option(
  "--generate-latex", dest="generateLatex", type="string",
  help="Generate the Latex (*.tex) output file using provided script" \
    " (i.e. rst2latex)",
  default="rst2latex.py" )
  
clp.add_option(
  "--generate-pdf", dest="generatePDF", type="string",
  help="Generate the PDF output file from the latex file using provided" \
    " script (i.e. pdflatex)",
  default="pdflatex" )

(options, args) = clp.parse_args(sys.argv)

#
# D) Get information for everything
#

if options.projectBaseDir:
  projectBaseDir = options.projectBaseDir
else:
  projectBaseDir = \
    os.path.abspath(os.path.join(tribitsBaseDir, "../.."))
#print "projectBaseDir =", projectBaseDir

projectName = getProjectName(projectBaseDir)
#print "projectName =", projectName

if options.outputFileBase:
  outputFileBase = options.outputFileBase
else:
  outputFileBase = os.path.join(projectBaseDir,
    projectName+"BuildQuickRef")
#print "outputFileBase =", outputFileBase

#
# E) Read in standard body and make substitution
#

tribitsBuildQuickRefBodyFile = \
  os.path.join(thisFileRealAbsBasePath, "TribitsBuildQuickRefBody.rst")
#print "tribitsBuildQuickRefBodyFile =", tribitsBuildQuickRefBodyFile

tribitsBuildQuickRefBodyStr = \
  readStrFromFile(tribitsBuildQuickRefBodyFile)

substitutedTribitsBuildQuickRefBodyStr = \
  tribitsBuildQuickRefBodyStr.replace("<Project>", projectName)

#
# F) Generate the output files
#

filesToClean = []

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

outputFileBaseName = os.path.basename(outputFileBase)

outputRstFile = outputFileBase+".rst"
print "Writing rst file ..."
openWriteFilePermissions(outputRstFile)
writeStrToFile(outputRstFile, projectBuildQuickRefStr)
setGeneratedFilePermissions(outputRstFile)

if options.generateHtml:
  print "Generating "+outputFileBaseName+".html ..."
  outputHtmlFile = outputFileBase+".html"
  generateFile(outputHtmlFile,
    options.generateHtml+" "+outputRstFile+" "+outputHtmlFile)

if options.generateLatex:
  print "Generating "+outputFileBaseName+".tex ..."
  outputLatexFile = outputFileBase+".tex"
  runSysCmnd(options.generateLatex+" "+outputRstFile+" "+outputLatexFile)
  if options.generatePDF:
    print "Generating "+outputFileBaseName+".pdf ..."
    outputPdfFile = outputFileBase+".pdf"
    outputPdfFileLog = outputLatexFile+".log"
    generateFile(outputPdfFile,
      options.generatePDF+" "+outputLatexFile,
      outFile=outputPdfFileLog,
      workingDir=projectBaseDir,
      runTwice=True)
    filesToClean.append(outputPdfFileLog)

#
# G) Clean the intermediate files
#

print "Cleaning intermediate files ..."

filesToClean.extend(
  [
    outputFileBase+".aux",
    outputFileBase+".log",
    outputFileBase+".out",
    outputFileBase+".tex",
    outputFileBase+".toc",
    ]
  )

for tempFile in filesToClean:
  if os.path.exists(tempFile):
    runSysCmnd("rm "+tempFile)
