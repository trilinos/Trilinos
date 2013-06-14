#!/usr/bin/env python

usageHelp = r"""create-project-user-quickstart.py [options]

Generate a TriBITS project-specific user quickstart file.

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
thisFileRealAbsBasePath = os.path.dirname(os.path.abspath(os.path.realpath(thisFilePath)))
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

def generateFile(filePath, generateCmnd, outFile=None, workingDir=""):
  openWriteFilePermissions(filePath)
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
  help="Project-specific user quickstart file base name." \
    "  If not set, then set to <projectBaseDir>/<Project>UserQuickstart"
  )
  
clp.add_option(
  "--generate-html", dest="generateHtml", type="string",
  help="Generate the HTML output file using provided script (i.e. rst2html)",
  default="rst2html" )
  
clp.add_option(
  "--generate-latex", dest="generateLatex", type="string",
  help="Generate the Latex (*.tex) output file using provided script" \
    " (i.e. rst2latex)",
  default="rst2latex" )
  
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
  outputFileBase = os.path.join(projectBaseDir, projectName+"UserQuickstart")
#print "outputFileBase =", outputFileBase

#
# E) Read in standard body and make substitution
#

tribitsUserQuickstartBodyFile = \
  os.path.join(thisFileRealAbsBasePath, "TribitsUserQuickstartBody.rst")
#print "tribitsUserQuickstartBodyFile =", tribitsUserQuickstartBodyFile

tribitsUserQuickstartBodyStr = \
  readStrFromFile(tribitsUserQuickstartBodyFile)

substitutedTribitsUserQuickstartBodyStr = \
  tribitsUserQuickstartBodyStr.replace("<Project>", projectName)

#
# F) Generate the output files
#

projectUserQuickstartTemplateFile = \
  os.path.join(projectBaseDir, "cmake", \
    projectName+"UserQuickstartTemplate.rst")
#print "projectUserQuickstartTemplateFile =", projectUserQuickstartTemplateFile

projectUserQuickstartTemplateStr = \
  readStrFromFile(projectUserQuickstartTemplateFile)

projectUserQuickstartStr = \
  projectUserQuickstartTemplateStr \
  + "\n\n" \
  + substitutedTribitsUserQuickstartBodyStr

outputFileBaseName = os.path.basename(outputFileBase)

outputRstFile = outputFileBase+".rst"
print "Writing rst file ..."
openWriteFilePermissions(outputRstFile)
writeStrToFile(outputRstFile, projectUserQuickstartStr)
setGeneratedFilePermissions(outputRstFile)

if options.generateHtml:
  print "Generating "+outputFileBaseName+".html ..."
  outputHtmlFile = outputFileBase+".html"
  generateFile(outputHtmlFile,
    options.generateHtml+" "+outputRstFile+" "+outputHtmlFile)

if options.generateLatex:
  print "Generating "+outputFileBaseName+".tex ..."
  outputLatexFile = outputFileBase+".tex"
  generateFile(outputLatexFile, \
    options.generateLatex+" "+outputRstFile+" "+outputLatexFile)
  if options.generatePDF:
    print "Generating "+outputFileBaseName+".pdf ..."
    outputPdfFile = outputFileBase+".pdf"
    generateFile(outputPdfFile,
      options.generatePDF+" "+outputLatexFile,
      outFile=outputLatexFile+".log",
      workingDir=projectBaseDir)
