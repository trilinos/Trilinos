#!/bin/env python3

usageHelp = r"""

Replace a given token string with another string in a file (but only touch the
file if there were changes).

This will only match complete tokens that match the regex:

  ([^A-Za-z0-9_])|^)<token-to-replace)[^A-Za-z0-9_]
"""


def getCmndLineOptions():
  from argparse import ArgumentParser, RawDescriptionHelpFormatter

  clp = ArgumentParser(description=usageHelp,
    formatter_class=RawDescriptionHelpFormatter)

  clp.add_argument(
    "-t", dest="tokenToReplace", required=True,
    help="Token to repalce" )

  clp.add_argument(
    "-r", dest="replacementString", required=True,
    help="Replacement string" )

  clp.add_argument(
    "-f", dest="inputFile", required=True,
    help="Input file (and also output if -o <file> not specified)" )

  clp.add_argument(
    "-o", dest="outputFile", default="",
    help="Output file (optional)" )

  options = clp.parse_args(sys.argv[1:])

  if options.outputFile == "":
    options.outputFile = options.inputFile

  return options


#
#  Main()
#

if __name__ == '__main__':

  import sys, re

  inOptions = getCmndLineOptions()

  beginLineTokenPattern = re.compile(
    fr'^{inOptions.tokenToReplace}([^A-Za-z0-9_])' )
  midLineTokenPattern = re.compile(
    fr'([^A-Za-z0-9_]){inOptions.tokenToReplace}([^A-Za-z0-9_])' )

  with open(inOptions.inputFile, 'r') as file:
    lines = file.readlines()

  fileWasChanged = False
  newLines = []
  for line in lines:
    newLine = beginLineTokenPattern.sub(
      inOptions.replacementString + r'\1',
      line)
    newLine = midLineTokenPattern.sub(
      r'\1' + inOptions.replacementString + r'\2',
      newLine)
    if newLine != line:
      fileWasChanged = True
    newLines.append(newLine)

  if (fileWasChanged or inOptions.outputFile != inOptions.inputFile):
    with open(inOptions.outputFile, 'w') as file:
      file.writelines(newLines)
