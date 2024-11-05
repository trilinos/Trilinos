#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import re


# Define regexes for matches

commandCallRegex = \
  re.compile(r'[a-zA-Z_][a-zA-Z0-9_]*\(')
commandCallWithSpacesRegex = \
  re.compile(r'[a-zA-Z_][a-zA-Z0-9_]*[\s]*\(')
macroDefRegex = \
  re.compile(r'[^a-zA-Z_][mM][aA][cC][rR][oO]\([\s]*[a-zA-Z_][a-zA-Z0-9_]*')
functionDefRegex = \
  re.compile(r'[^a-zA-Z_][fF][uU][nN][cC][tT][iI][oO][nN]\([\s]*[a-zA-Z_][a-zA-Z0-9_]*')

specialCommandNamesList = [
  r'[sS][eE][tT]',
  r'[iI][fF]',
  r'[eE][lL][sS][eE][iI][fF]',
  r'[fF][oO][rR][eE][aA][cC][hH]',
  r'[wW][hH][iI][lL][eE]',
  ]
specialCommandRegexList = []
for specialCommandName in specialCommandNamesList:
  specialCommandRegexList.append(
    re.compile(specialCommandName+r'[\s]*\(') )

andLogicalRegex = re.compile(r'AND[\s]*\(', re.DOTALL)
orLogicalRegex = re.compile(r'OR[\s]*\(', re.DOTALL)
notLogicalRegex = re.compile(r'NOT[\s]*\(', re.DOTALL)


# Lower case CMake command calls and macro and function names.
#
def makeCmndsLowerCaseInCMakeStr(cmakeCodeStrIn):
  cmakeCodeStrOut = cmakeCodeStrIn
  cmakeCodeStrOut = makeRegexMatchesLowerCaseInCMakeStr(cmakeCodeStrOut,
    commandCallWithSpacesRegex)
  cmakeCodeStrOut = makeRegexMatchesLowerCaseInCMakeStr(cmakeCodeStrOut,
    macroDefRegex)
  cmakeCodeStrOut = makeRegexMatchesLowerCaseInCMakeStr(cmakeCodeStrOut,
    functionDefRegex)
  return cmakeCodeStrOut


# Make regex pattern matches lower cases.
#
def makeRegexMatchesLowerCaseInCMakeStr(cmakeCodeStrIn, regex):
  cmakeCodeStrOut = ""
  cmakeCodeStrInLen = len(cmakeCodeStrIn)
  cmakeCodeStrStartSearchIdx = 0
  while cmakeCodeStrStartSearchIdx < cmakeCodeStrInLen:
    m = regex.search(cmakeCodeStrIn[cmakeCodeStrStartSearchIdx:])
    if m:
      #print("")
      #print("m.group() = '"+str(m.group())+"'")
      # Get the substr for the match
      cmakeCodeStrMatchStartIdx = cmakeCodeStrStartSearchIdx + m.start() 
      cmakeCodeStrMatchEndIdx = cmakeCodeStrStartSearchIdx + m.end() 
      commandCallMatchStr = \
        cmakeCodeStrIn[cmakeCodeStrMatchStartIdx:cmakeCodeStrMatchEndIdx]
      assert m.group() == commandCallMatchStr
      # Lower case the matching command call and update the output string
      # since last match
      cmakeCodeStrOut += \
        cmakeCodeStrIn[cmakeCodeStrStartSearchIdx:cmakeCodeStrMatchStartIdx]
      cmakeCodeStrOut += conditionalLowerCaseCMakeGroup(commandCallMatchStr)
      # Update indexes for next search
      cmakeCodeStrStartSearchIdx += m.end()
    else:
      # No more matches so copy the rest of the string
      cmakeCodeStrOut += \
        cmakeCodeStrIn[cmakeCodeStrStartSearchIdx:]
      break

  return cmakeCodeStrOut


def conditionalLowerCaseCMakeGroup(cmakeMatchGroup):
  cmakeMatchGroupNewCase = None
  #print("")
  #print("cmakeMatchGroup = '"+cmakeMatchGroup+"'")
  if isConditionalCMakeGroup(cmakeMatchGroup):
    #print("Is a conditional expression")
    cmakeMatchGroupNewCase = cmakeMatchGroup
  elif isSpecialCommandNameMatch(cmakeMatchGroup):
    #print("Is special command name!")
    cmakeMatchGroupNewCase = cmakeMatchGroup.lower()
  elif isMacroOrFunctionGroup(cmakeMatchGroup):
    cmakeMatchGroupNewCase = cmakeMatchGroup.lower()
  elif commandCallRegex.match(cmakeMatchGroup):
    #print("Is standard command call!")
    cmakeMatchGroupNewCase = cmakeMatchGroup.lower()
  else:
    #print("Does not match known cases so don't lower case!")
    cmakeMatchGroupNewCase = cmakeMatchGroup
  #print("cmakeMatchGroupNewCase = '"+cmakeMatchGroupNewCase+"'")
  return cmakeMatchGroupNewCase


def isConditionalCMakeGroup(cmakeMatchGroup):
  if andLogicalRegex.match(cmakeMatchGroup) \
    or orLogicalRegex.match(cmakeMatchGroup) \
    or notLogicalRegex.match(cmakeMatchGroup) \
    :
    return True
  return False


def isSpecialCommandNameMatch(cmakeMatchGroup):
  for specalCommandRegex in specialCommandRegexList:
    #print("specalCommandRegex.pattern = '"+specalCommandRegex.pattern+"'")
    m = specalCommandRegex.match(cmakeMatchGroup)
    if m:
      #print("MATCH!")
      return True
  return False


def isMacroOrFunctionGroup(cmakeMatchGroup):
  if macroDefRegex.match(cmakeMatchGroup):
    return True
  if functionDefRegex.match(cmakeMatchGroup):
    return True
  return False


#
# Helper functions for main()
#


usageHelp = r"""

Convert cmake command calls to lower case and macro and function names in
definitions to lower case in a given file.

The replacements are somewhat conservative in order not to make too many
replacements in non-CMake code.  That is, in general, any text of the form:

  <identifier>(

will result in '<identifier>' being lower-cased.  However, this does not apply
for logical operators 'AND, 'OR', and 'NOT'.

Also, certain common CMake commands that often have spaces between
'<identifier>' and '(' like:

  SET (...)
  IF (...)
  ELSEIF (...)
  FOREACH (...)
  WHILE (...)

will be lower-cased as well.  (These are less likely to match in regular
text.)
"""


def getCmndLineOptions():
  from argparse import ArgumentParser, RawDescriptionHelpFormatter
  clp = ArgumentParser(description=usageHelp,
    formatter_class=RawDescriptionHelpFormatter)
  clp.add_argument("file",
    help="The file to lower-case the CMake commands within." )
  options = clp.parse_args(sys.argv[1:])
  return options


def validateCmndLineOptions(inOptions):
  if not os.path.exists(inOptions.file):
    print("Error, the file '"+inOptions.file+"' does not exist!")
    sys.exit(1)


def lowerCaseCMakeCodeInFile(inOptions):
  with open(inOptions.file, 'r') as fileHandle:
    cmakeStrOrig = fileHandle.read()
  cmakeStrLowerCased = makeCmndsLowerCaseInCMakeStr(cmakeStrOrig)
  with open(inOptions.file, 'w') as fileHandle:
    fileHandle.write(cmakeStrLowerCased)

#
# main()
#


if __name__ == '__main__':
  inOptions = getCmndLineOptions()
  validateCmndLineOptions(inOptions)
  lowerCaseCMakeCodeInFile(inOptions)
