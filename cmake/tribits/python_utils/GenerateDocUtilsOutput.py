import sys
import os
import stat
import subprocess

#
# A) Set up basic paths and import modules
#

from GeneralScriptSupport import *


#
# B) Define some helper functions
#


def openWriteFilePermissions(filePath):
  if os.path.exists(filePath):
    os.chmod(filePath, stat.S_IREAD | stat.S_IWRITE \
      | stat.S_IRGRP | stat.S_IWGRP)


def setGeneratedFilePermissions(filePath):
  os.chmod(filePath, stat.S_IREAD | stat.S_IRGRP | stat.S_IROTH)


def generateFile(filePath, generateCmnd, outFile=None, workingDir="",
  runTwice=False, echoCmnds=False \
  ):
  openWriteFilePermissions(filePath)
  runSysCmnd(generateCmnd, outFile=outFile, workingDir=workingDir,
    echoCmndForDebugging=echoCmnds)
  if runTwice:
    runSysCmnd(generateCmnd, outFile=outFile, workingDir=workingDir,
       echoCmndForDebugging=echoCmnds)
  setGeneratedFilePermissions(filePath)


def addCmndLineOptions(clp):

  # Find the right default for the current system
  rst2html = "rst2html"
  rst2latex = "rst2latex"
  rst2htmlWhich = getCmndOutput("which rst2html", True, False)
  if rst2htmlWhich == "" or re.match(".+no rst2html.+", rst2htmlWhich):
    rst2html = rst2html+".py"
    rst2latex = rst2latex+".py"

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
      " [Default '"+rst2html+"']",
    default=rst2html )
    
  clp.add_option(
    "--generate-latex", dest="generateLatex", type="string",
    help="Generate the Latex (*.tex) output file using provided script" \
      " (i.e. rst2latex) [Default '"+rst2latex+"']",
    default=rst2latex )
    
  clp.add_option(
    "--generate-latex-options", dest="generateLatexOptions", type="string",
    help="Options to pass to the generate latex command",
    default="" )
    
  clp.add_option(
    "--generate-pdf", dest="generatePDF", type="string",
    help="Generate the PDF output file from the latex file using provided" \
      " script (i.e. pdflatex) [Default 'pdflatex']",
    default="pdflatex" )

  clp.add_option(
    "--clean-temp-files", dest="cleanTempFiles", action="store_true",
    help="Clean temporary files used in generation. [default]" )
  clp.add_option(
    "--no-clean-temp-files", dest="cleanTempFiles", action="store_false",
    help="Do not delete temporary files.",
    default=True )

  clp.add_option(
    "--echo-cmnds", dest="echoCmnds", action="store_true",
    help="Echo the commands being run to STDOUT." )
  clp.add_option(
    "--no-echo-cmnds", dest="echoCmnds", action="store_false",
    help="Do not echo the commands being run.",
    default=False )


def generateDocutilsOuputFiles(options):

  filesToClean = []
  
  # Base name including path (must just be relative)
  outputFileBase = options.fileBase
  # The path of the rst file:
  rstFile = outputFileBase+".rst"
  # Just the base name
  outputFileBaseName = os.path.basename(outputFileBase)
  
  if options.generateHtml:
    print("Generating " + outputFileBaseName + ".html ...")
    outputHtmlFile = outputFileBase+".html"
    generateFile(outputHtmlFile,
      options.generateHtml+" "+rstFile+" "+outputHtmlFile, echoCmnds=options.echoCmnds)
  
  if options.generateLatex:
    print("Generating " + outputFileBaseName + ".tex ...")
    outputLatexFile = outputFileBase+".tex"
    runSysCmnd(options.generateLatex+" "+options.generateLatexOptions+ \
       " "+rstFile+" "+outputLatexFile, echoCmndForDebugging=options.echoCmnds)
    if options.generatePDF:
      print("Generating " + outputFileBaseName + ".pdf ...")
      outputPdfFile = outputFileBase+".pdf"
      outputPdfFileLog = outputLatexFile+".log"
      generateFile(outputPdfFile,
        options.generatePDF+" -halt-on-error "+outputLatexFile,
        outFile=outputPdfFileLog,
        runTwice=True, echoCmnds=options.echoCmnds)
      filesToClean.append(outputPdfFileLog)
  
  #
  # Clean the intermediate files
  #
  
  if options.cleanTempFiles:

    print("Cleaning intermediate files ...")
    
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

  else:

    print("Keeping temp files ...")
