#!/usr/bin/env python

# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
# @HEADER

import sys
import os
import subprocess
import re
import pprint
import glob


pp = pprint.PrettyPrinter(indent=4)


#
# Script documentation
#


usageHelp = r"""extract_rst_cmake_doc.py [OPTIONS]

This script implements a system to extract RST-formatted (reStructuredText)
comment blocks out of *.cmake files and populate them in a set of template
*.rst files.

An example of usage would be::

$ cd tribits/doc/developers_guide
$ extract_rst_cmake_doc.py \
  --extract-from=../../package_arch/,../..//utils/ \
  --rst-file-pairs=TribitsDetailedMacroFunctionDocsTemplate.rst:TribitsDetailedMacroFunctionDocs.rst

The --extract-from argument will read every *.cmake file in the
directory given if the directory is specified with '/' at the end.  We use
--file-extensions (default is '.cmake') to provide the set of file patterns.

The --rst-file-pairs argument gives pairs of files
<name>Template.rst:<name>.rst that will subsutute the placeholders for macros
and function definitions from the file <name>Template.rst and will produce the
file <name>.rst.

The way it works is that the *.cmake files specified by
--extract-from must contain comment blocks of the form:

{{{
#
# @MACRO: SOME_MACRO_NAME()
#
# [rst formatted macro doc]
#
MACRO(SOME_MACRO_NAME ...
...

#
# @FUNCTION: SOME_FUNC_NAME()
#
# [rst formatted func doc]
#
FUNCTION(SOME_FUNC_NAME ...
...
}}}

The rules for these comment blocks are:

* The comment blocks '#' must begin on the first column.

* The special comment block identifier '@MACRO: <name>()' must start on the
  3rd column, i.e. '# @MACRO: <name>()'.  Otherwise, the comment block will be
  ignored.

* Text in the comment blocks must begin after the 2nd column.  That is,
  '#somthing' is illegal and will result in an error.

* The comment blocks must terminate at the MACRO or FUNCTION definition and
  the <name> of the macro or function listed in the marker '@MACRO: <name>()'
  or '@MACRO: <name>()' must match what is in the MACRO(<name> ...) or
  FUNCTION(<name> ...)  definition.

* Whitespce is allowed in the block marker between '@MACRO' and ':' and after
  the definition 'MACRO(' and 'FUNTION('.  However, the name must be on the
  same line.

These RST comment blocks are then inserted into RST template files
<name>Template.rst producing output files <name>.rst listed in the
--rst-file-pairs argument.

The format of these RST template files is:

{{{

Detailed Macro and Function Documentation
=========================================

@FUNCTION: SOME_FUNC_NAME2() -

@MACRO: SOME_MACRO_NAME1() -

}}}

The line of the form '@FUNCTION: SOME_FUNC_NAME2() -' are flagged as
replacement blocks that are substitted with the extracted RST blocks of the
same name and type.  For example, the above RST template file would be
substituted to look like:

{{{

Detailed Macro and Function Documentation
=========================================

SOME_FUNC_NAME2()
-----------------

[rst formatted func text]


SOME_MACRO_NAME1()
------------------

[rst formatted macro text]


}}}

The format of these substitution lines is:

  @<blockType>: <blockName> <sepChar>

where:

* <blockType> is either 'MACRO' or 'FUNCTION'.
* <blockName> is the function or macro name (including the parath).
* <sepChar> is a single char used to build the RST section underline.

The rules for this format are:

* The '@' in '@<blockType> must being on the first column.  Otherwise the
  substitution block is ignored.

* <blockName> must match the name of a previously read in RST documentation
  block.  If it is not, then processing ends with an error.

* The <sepChar> must be present and must be a single char, typically '-',
  '+', "=", etc.  This should be picked to be consistent with the RST document
  this is being put into.

* At least one space must separate the fields '<blockName> <sepChar>' or an
  error will occur.

See the unit tests in extract_rst_cmake_doc_UnitTest.py for examples of
behavior w.r.t. different errors and formatting mistakes.
"""

#
# Helper functions
#


# Strip off the "()" for the entity name
def getBaseEntityName(entityName):
  parenthIdx = entityName.find("(")
  if parenthIdx > 0:
    return entityName[0:parenthIdx]
  return entityName


# Get the line type and entity name for CMake code
def getLineEntityTypeAndName(line, fileNameAndLinePrefix, currentRstDocBlockType,
  inRstDocBlockBaseName \
  ):
  lineSplit = line.split("(")
  if len(lineSplit) == 1:
    raise Exception(fileNameAndLinePrefix+ \
      "error: expecting "+currentRstDocBlockType+"("+inRstDocBlockBaseName+" ...)"+\
      " on this line.  RST comment block must terminate in the stated entity!")
  lineEntityType = lineSplit[0].strip()
  lineEntityName = lineSplit[1].strip().split(" ")[0].strip()
  if lineEntityName[-1] == ")":
    lineEntityName = lineEntityName[0:-1]
  return (lineEntityType.upper(), lineEntityName)


# Format file name and line number
def getFileNameLineNum(fileName, lineNum):
  return fileName+":"+str(lineNum)


# Format file name and line number as prefix
def getFileNameLineNumPrefix(fileName, lineNum):
  return getFileNameLineNum(fileName, lineNum)+": "


# Extract a set of RST doc blocks from the given text block (read from a file).
#
def extractRstDocBlocksFromText(rawText, rstBlockTypes, fileName, traceExtraction=False \
  ):

  rstDocBlocks = {}

  inCommentBlock = False

  inRstDocBlock = ""
  currentRstDocBlockType = ""
  currentRstBlockBody = ""

  line_i = 0
  for line in rawText.splitlines():
    #print("\nline = '" + line + "'")

    line_i += 1

    fileNameAndLinePrefix = fileName+":"+str(line_i)+": "

    # Determine if this is a comment line
    if len(line) == 0:
      inCommentBlock = False
    elif line[0] == '#':
      inCommentBlock = True
    else:
      inCommentBlock = False

    justFoundBlockBeginning = False
    if inCommentBlock and not inRstDocBlock:
      # Look for beginning of expected RST doc block type
      for blockType in rstBlockTypes:
        blockBeginning = "# @"+blockType
        lineBlockBeginning = line[0:len(blockBeginning)]
        #print("blockBeginning = '" + blockBeginning + "'")
        #print("lineBlockBeginning = '" + lineBlockBeginning + "'")
        if lineBlockBeginning == blockBeginning:
          #print("Matches block beginning " + blockType + "!")
          currentRstDocBlockType = blockType
          splitOnColon = line.split(":")
          if len(splitOnColon) != 2:
            raise Exception(fileNameAndLinePrefix+ \
              "error: '"+line+"' is missing the colon ':' separator!")
          inRstDocBlock = splitOnColon[1].strip()
          #print("inRstDocBlock = '" + inRstDocBlock + "'")
          #print("currentRstDocBlockType = '" + currentRstDocBlockType + "'")
          inRstDocBlockFileNameLineNum = getFileNameLineNum(fileName, line_i)
          if traceExtraction:
            print("Extracting '" + blockType + "' block '" + inRstDocBlock + "'"
                  + " from " + getFileNameLineNumPrefix(fileName, line_i))
          currentRstBlockBody = ""
          justFoundBlockBeginning = True
          break
   
    if justFoundBlockBeginning:
      # Just found the beginning of the RST block so move to next line!
      continue

    if inCommentBlock and inRstDocBlock:
      # Strip off the beginning "# " and add to body
      if len(line) > 1 and line[1] != " ":
        raise Exception( \
          getFileNameLineNumPrefix(fileName, line_i)+ \
          "error: Comment blocks must have at least one space after '#'"+\
          " in line = '"+line+"'")
      rstBlockBodyLine = line[2:]
      #print("rstBlockBodyLine = '" + rstBlockBodyLine + "'")
      currentRstBlockBody += (rstBlockBodyLine + "\n")
      #print("currentRstBlockBody = {\n" + currentRstBlockBody + "}")

    if not inCommentBlock and inRstDocBlock:
      # This is the end of the RST block
      fileNameAndLinePrefix = getFileNameLineNumPrefix(fileName, line_i)
      inRstDocBlockBaseName = getBaseEntityName(inRstDocBlock)
      # First, verify that the comment block matches the entity
      (lineEntityType, lineEntityName) = getLineEntityTypeAndName(
         line, fileNameAndLinePrefix, currentRstDocBlockType, inRstDocBlockBaseName)
      if lineEntityType != currentRstDocBlockType:
        raise Exception(fileNameAndLinePrefix+ \
          "error: expecting '"+currentRstDocBlockType+\
          "' but got type '"+lineEntityType+"''")
      if lineEntityName != inRstDocBlockBaseName:
        raise Exception(fileNameAndLinePrefix+ \
          "error: expecting "+lineEntityType+" "+inRstDocBlock+\
          " but got wrong "+lineEntityType+" name '"+lineEntityName+"'")
      # ToDo: Assert that blockName is the same!
      # Next, story the RST block
      rstDocBlocks.update( {
        inRstDocBlock : {
          "type" : currentRstDocBlockType,
          "body" : currentRstBlockBody,
          "fileNameLineNum" : inRstDocBlockFileNameLineNum,
          }
        } )
      # Finally, rest to look for the next RST block
      inRstDocBlock = ""
      currentRstDocBlockType = ""
      #print("rstDocBlocks:\n" + str(rstDocBlocks))
      
      # ToDo: Check that this line starts with currentRstDocBlockType or
      # throw exception!

  # end for line

  return rstDocBlocks


#
# Get an RST block from dict or throw except if it does not exist
#
def getRstDocBlock(rstDocBlocks, rstBlockName, fileNameAndLinePrefix, line):
  rstDocBlock = rstDocBlocks.get(rstBlockName, "")
  if rstDocBlock == "":
    raise Exception(fileNameAndLinePrefix+"Error, block name"+\
      " '"+rstBlockName+"' does not exist in rstDocBlocks for line = '"+line+"'")
  return rstDocBlock


#
# Get RST section separator of length lineLen
#
def getRstSectStr(sepChar, lineLen):
  return sepChar * lineLen


def removeEmtpyElements(arrayIn):
  arrayOut = []
  for ele in arrayIn:
    if ele != "":
      arrayOut.append(ele)
  return arrayOut


# Replace references to RST blocks with sections
#
def replaceWithRstDocBlocksInText(textToReplace, rstBlockTypes, rstDocBlocks, fileName, \
  traceReplacements=False, fileNamePathBaseDir="", includeFileNameLineNum=False \
  ):

  replacedText = ""

  textToReplaceLines = textToReplace.splitlines()
  numTextToReplaceLines = len(textToReplaceLines)

  if not fileNamePathBaseDir == "":
    fileNamePathBaseDirAndSepLen = len(fileNamePathBaseDir)
    if not fileNamePathBaseDir.endswith("/"): fileNamePathBaseDirAndSepLen += 1
  else:
    fileNamePathBaseDirAndSepLen = 0

  line_i = 0

  for line in textToReplaceLines:

    line_i += 1

    #print("\nline_i = " + str(line_i))
    #print("line = '" + line + "'")

    # Look for block of text to replace 
    replacedRstBlock = False
    for blockType in rstBlockTypes:
      blockBeginning = "@"+blockType
      lineBlockBeginning = line[0:len(blockBeginning)]
      #print("blockBeginning = '" + blockBeginning + "'")
      #print("lineBlockBeginning = '" + lineBlockBeginning + "'")
      if lineBlockBeginning == blockBeginning:
        #print("Matches block beginning " + blockType + "!")
        lineSplitArray = line.split(":")
        fileNameAndLinePrefix = fileName+":"+str(line_i)+": "
        #print("lineSplitArray = " + str(lineSplitArray))
        if len(lineSplitArray)==1:
          # Missing the colon so no replacement!
          break
        rstBlockNameAndSecCharArray = removeEmtpyElements(
          lineSplitArray[1].strip().split(" "))
        #print("rstBlockNameAndSecCharArray = " + str(rstBlockNameAndSecCharArray))
        if len(rstBlockNameAndSecCharArray) != 2:
          raise Exception(fileNameAndLinePrefix+"Error, line does not match format"+\
            " '@"+blockType+": <blockName> <sepChar>' for line = '"+line+"'")
        rstBlockName = rstBlockNameAndSecCharArray[0].strip()
        rstBlockSecChar = rstBlockNameAndSecCharArray[1].strip()
        #print("rstBlockName = '" + rstBlockName + "'")
        #print("rstBlockSecChar = '" + rstBlockSecChar + "'")
        if len(rstBlockSecChar) != 1:
          raise Exception(fileNameAndLinePrefix+\
            "Error, separation char must be on char in line = '"+line+"'")
        if traceReplacements:
            print("Replacing '" + blockType + "' block '" + rstBlockName + "'" +
                  " in " + getFileNameLineNumPrefix(fileName, line_i))
        rstDocBlock = getRstDocBlock(rstDocBlocks, rstBlockName,
          fileNameAndLinePrefix, line)
        if rstDocBlock.get("type") != blockType:
          raise Exception(fileNameAndLinePrefix+"Error, looked up block type"+\
            " '@"+rstDocBlock.get("type")+"' does not match block type for line = '"+line+"'")
        rstSecStr = getRstSectStr(rstBlockSecChar, len(rstBlockName))
        # Add the replaced block section
        replacedText += (rstBlockName+"\n")
        replacedText += (rstSecStr+"\n")
        replacedText += rstDocBlock.get("body")
        if includeFileNameLineNum:
          fileNameLineNum = rstDocBlock.get("fileNameLineNum")
          if fileNamePathBaseDirAndSepLen > 0:
            relFileNameLineNum = fileNameLineNum[fileNamePathBaseDirAndSepLen:]
          else:
            relFileNameLineNum = fileNameLineNum
          replacedText += "\nIn: "+relFileNameLineNum + "\n\n"
        replacedRstBlock = True
        break

    if not replacedRstBlock and not (line_i == numTextToReplaceLines and line==""):
      replacedText += (line+"\n")

    #print("replacedText = {\n" + replacedText + "}")

  # end for line

  return replacedText


# Get a list of files give list of files, dirs, and extensions
def getExtractFilesList(inOptions):
  filesList = []
  for fileOrDir in inOptions.extractFrom.split(","):
    if fileOrDir[-1] == "/":
      # This is a directory so grab all files with the matching extensions
      for fileExtension in inOptions.fileExtensions.split(","):
        globStr = fileOrDir+"*"+fileExtension
        globedFiles = glob.glob(globStr)
        filesList.extend(globedFiles)
    else:
      filesList.append(fileOrDir)
  return filesList


# Get the list of template and replacement files
def getRstFilesList(inOptions):
  rstFilesList = []
  if inOptions.rstFilePairs:
    for templateFileAndFile in inOptions.rstFilePairs.split(","):
      rstFilesList.append(templateFileAndFile.split(":"))
  return rstFilesList


# Extract RST blocks out of a set of files
def extractRstDocBlocksFromFileList(fileList, rstBlockTypes, traceExtraction=False):
  rstDocBlocks = {}
  for fileName in fileList:
    rawFileStr =   open(fileName, 'r').read()
    rstDocBlocks.update(
      extractRstDocBlocksFromText(rawFileStr, rstBlockTypes, fileName, traceExtraction)
      )
  return rstDocBlocks


# Make substitutions in a template file and produce a target file
def replaceWithRstDocBlocksInTemplateFileList(rstFilesList, rstBlockTypes,
  rstDocBlocks, traceReplacements=False, fileNamePathBaseDir="",
  includeFileNameLineNum=False \
  ):
  for (templateFileName, fileName) in rstFilesList:
    templateFileStr =   open(templateFileName, 'r').read()
    fileStr = replaceWithRstDocBlocksInText(templateFileStr, rstBlockTypes,
      rstDocBlocks, templateFileName, traceReplacements=traceReplacements,
      fileNamePathBaseDir=fileNamePathBaseDir,
      includeFileNameLineNum=includeFileNameLineNum )
    open(fileName, 'w').write(fileStr)
 
  
# Run a command and synchronize the output
def runCmnd(options, cmnd):
  if options.debug:
    print("*** Running command: " + cmnd)
  if options.noOpt:
    print(cmnd)
  else:
    child = subprocess.Popen(cmnd, stdout=subprocess.PIPE).stdout
    output = child.read()
    sys.stdout.flush()
    print(output)
    sys.stdout.flush()


#
# Script body
#


if __name__ == '__main__':

  #
  # A) Get command-line options
  #
  
  from optparse import OptionParser
  
  clp = OptionParser(usage=usageHelp)

  clp.add_option(
    "--extract-from", dest="extractFrom", type="string",
    default="",
    help="List (comma separated) of directories (ending with '/') and files that the RST"+\
    " comment blocks will be extracted from.  Directories ending with '/' will result in"+\
    " files being found with extensions listed in --file-extensions.  Otherwise, the given"+\
    " files are searched.")

  clp.add_option(
    "--file-extensions", dest="fileExtensions", type="string",
    default=".cmake",
    help="List (comma separated) of file extensions (e.g. '.cmake').  The default"+\
    "is '.cmake'")

  clp.add_option(
    "--rst-file-pairs", dest="rstFilePairs", type="string",
    default="",
    help="List (comma separated) files of the form '<name>Template.rst:<name>.rst,...'"+\
    " that gives sets of template files to have comment blocks inserted into.")

  clp.add_option(
    "--show-file-name-line-num", dest="showFileNameLineNum", action="store_true",
    help="Include the relative file name path and line number for each documentation block.")
  clp.add_option(
    "--no-show-file-name-line-num", dest="showFileNameLineNum", action="store_false", default=False,
    help="Include the relative file name path and line number for each documentation block. [default]" )

  clp.add_option(
    "--file-name-path-base-dir", dest="fileNamePathBaseDir", type="string",
    default="",
    help="Base path stripped off of file names reported when --show-file-name-line-num is set."+\
      "  NOTE: This path should be relative to the paths in --extract-from and may be relative paths." )

  clp.add_option(
    "--dump-rst-blocks", dest="dumpRstBlocks", action="store_true",
    help="Print out the RST documentation blocks that are read in.")
  clp.add_option(
    "--no-dump-rst-blocks", dest="dumpRstBlocks", action="store_false", default=False,
    help="Print out the RST documentation blocks that are read in. [default]" )

  clp.add_option(
    "--do-trace", dest="doTrace", action="store_true",
    help="Print out trace of RST blocks read and replaced.")
  clp.add_option(
    "--no-do-trace", dest="doTrace", action="store_false", default=False,
    help="Do not print out trace of RST blocks read and replaced. [default]" )

  
  (options, args) = clp.parse_args()

  # Echo the command-line

  if options.doTrace:
    print("")
    print("**************************************************************************")
    print("Script: extract_rst_cmake_doc.py \\")
    print("  --extract-from='" + options.extractFrom + "' \\")
    print("  --file-extensions='" + options.fileExtensions + "' \\")
    print("  --rst-file-pairs='" + options.rstFilePairs + "' \\")
    if options.showFileNameLineNum:
      print("  --show-file-name-line-num \\")
    else:
      print("  --no-show-file-name-line-num \\")
    print("  --file-name-path-base-dir='" + options.fileNamePathBaseDir + "' \\")
    if options.dumpRstBlocks:
      print("  --dump-rst-blocks \\")
    else:
      print("  --no-dump-rst-blocks \\")
    if options.doTrace:
      print("  --do-trace \\")
    else:
      print("  --no-do-trace \\")

  # Assert the commandline

  if not options.extractFrom:
    print("\nError, --extract-from must be specified (see --help)!")
    sys.exit(1)

  if options.doTrace and not options.rstFilePairs:
    print("\nWarning: --rst-file-pairs is empty and no RST comment blocks will be set!")

  #
  # B) Read in all of the RST documentation blocks
  #

  rstBlockTypes = ["MACRO", "FUNCTION"]

  if options.doTrace:
    print("\nExtracting RST documentation blocks in --extract-from:")
  extractFromFilesList = getExtractFilesList(options)
  rstDocBlocks = extractRstDocBlocksFromFileList(
    extractFromFilesList, rstBlockTypes, options.doTrace)

  if options.dumpRstBlocks:
    print("Read in RST blocks:\n")
    pp.pprint(rstDocBlocks)
  
  #
  # C) Make the substitutions in all of the file pairs
  #

  if options.doTrace:
    print("\nReplacing RST documentation blocks in RST files in --rst-file-pairs:")
  rstFilesList = getRstFilesList(options)
  replaceWithRstDocBlocksInTemplateFileList(rstFilesList, rstBlockTypes,
    rstDocBlocks, traceReplacements=options.doTrace,
    fileNamePathBaseDir=options.fileNamePathBaseDir,
    includeFileNameLineNum=options.showFileNameLineNum )
