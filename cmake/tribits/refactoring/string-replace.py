#!/bin/env python3

usageHelp = r"""

Replace a given string with another string in a file (but only touch the file
if there were changes).

"""


def getCmndLineOptions():
  from argparse import ArgumentParser, RawDescriptionHelpFormatter

  clp = ArgumentParser(description=usageHelp,
    formatter_class=RawDescriptionHelpFormatter)

  clp.add_argument(
    "-s", dest="stringToReplace", required=True,
    help="String to repalce" )

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

  import sys

  inOptions = getCmndLineOptions()

  with open(inOptions.inputFile, 'r') as file:
    lines = file.readlines()

  fileWasChanged = False
  newLines = []
  for line in lines:
    newLine = line.replace(inOptions.stringToReplace, inOptions.replacementString)
    if newLine != line:
      fileWasChanged = True
    newLines.append(newLine)

  if (fileWasChanged or inOptions.outputFile != inOptions.inputFile):
    with open(inOptions.outputFile, 'w') as file:
      file.writelines(newLines)
