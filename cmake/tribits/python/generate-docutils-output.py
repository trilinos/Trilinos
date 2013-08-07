#!/usr/bin/env python

usageHelp = r"""generate-docutils-output.py --file-base=<FILEBASE> [options]

Create output files for DocUtils reStructuredText.

This simple script can be used to geenrate HTML, Latex, and PDF output for a
given <outputFileBase>.rst file.  What output is generated is controlled by
the options --generate-html, --generate-latex, and --generate-pdf.  This
requires that the Docutils commands be installed on the system as well as
latex (if generating PDF output).
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


if __name__ == '__main__':



  #
  # C) Read in the commandline options
  #
    
  clp = OptionParser(usage=usageHelp)
  
  clp.add_option(
    "--file-base", dest="fileBase", type="string",
    default="",
    help="Base name for the reStructuredText *.rst file.  This may include the" \
      " relative or absolute path up to but not including the '.rst' extension." \
      "  [Required]"
    )
    
  clp.add_option(
    "--generate-html", dest="generateHtml", type="string",
    help="Generate the HTML output file using provided script (i.e. rst2html)" \
      " [Default 'rst2html.py']",
    default="rst2html.py" )
    
  clp.add_option(
    "--generate-latex", dest="generateLatex", type="string",
    help="Generate the Latex (*.tex) output file using provided script" \
      " (i.e. rst2latex) [Default 'rst2latex.py']",
    default="rst2latex.py" )
    
  clp.add_option(
    "--generate-pdf", dest="generatePDF", type="string",
    help="Generate the PDF output file from the latex file using provided" \
      " script (i.e. pdflatex) [Default 'pdflatex']",
    default="pdflatex" )
  
  (options, args) = clp.parse_args(sys.argv)
  
  if not options.fileBase:
    print "Error, --file-base=<fileBase> must be specified!"
    sys.exit(1)
  
  #
  # D) Generate the output files
  #
  
  filesToClean = []
  
  # Base name including path (must just be relative)
  outputFileBase = options.fileBase
  # The path of the rst file:
  rstFile = outputFileBase+".rst"
  # Just the base name
  outputFileBaseName = os.path.basename(outputFileBase)
  
  if options.generateHtml:
    print "Generating "+outputFileBaseName+".html ..."
    outputHtmlFile = outputFileBase+".html"
    generateFile(outputHtmlFile,
      options.generateHtml+" "+rstFile+" "+outputHtmlFile)
  
  if options.generateLatex:
    print "Generating "+outputFileBaseName+".tex ..."
    outputLatexFile = outputFileBase+".tex"
    runSysCmnd(options.generateLatex+" "+rstFile+" "+outputLatexFile)
    if options.generatePDF:
      print "Generating "+outputFileBaseName+".pdf ..."
      outputPdfFile = outputFileBase+".pdf"
      outputPdfFileLog = outputLatexFile+".log"
      generateFile(outputPdfFile,
        options.generatePDF+" "+outputLatexFile,
        outFile=outputPdfFileLog,
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
  
