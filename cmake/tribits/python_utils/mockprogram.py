#!/usr/bin/env python

# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

#
# Usage: mockprogram.py [any arguments]
#
# Mock program that takes input arguments and produces stdout by reading from
# a file .mockprogram_inout.txt in the current directory or the file specified
# by the env var MOCKPROGRAM_INOUT_FILE_OVERRIDE (which can be in any
# directory).  This script is used to take the place of real commands during a
# test that involves calling commands on the commandline.
#
# The file .mockprogram_inout.txt (or pointed to by
# MOCKPROGRAM_INOUT_FILE_OVERRIDE) is of the form:
#
#   MOCK_PROGRAM_INPUT: <args_1>
#   MOCK_PROGRAM_RETURN: <rtn>
#   MOCK_PROGRAM_OUTPUT: <outline_1_line_1>
#   <outline_1_line_2>
#   ...
#   MOCK_PROGRAM_INPUT: <args_2>
#
# The program reads in the blocks starting at the time and removes the block
# from the file after it runs.  After all of the blocks are read in, if run
# again it will error out with error code 2.
#
# This program can be used, for example, to simulate git command.  For
# example, a couple of git commits might be simulated like:
#
#   MOCK_PROGRAM_INPUT: log -1
#   MOCK_PROGRAM_RETURN: 0
#   MOCK_PROGRAM_OUTPUT: This is the summary line
#
#   The is the body of the commit msg
#   MOCK_PROGRAM_INPUT: diff --name-only HEAD --not @{u}
#   MOCK_PROGRAM_RETURN: 0
#   MOCK_PROGRAM_OUTPUT: file_name_1.txt
#   file_name_2.txt
#   file_name_3.txt

#

import sys
import os

inputArgs = ' '.join(sys.argv[1:])
#print("inputArgs = '" + inputArgs + "'"

if os.environ.get("MOCKPROGRAM_INOUT_FILE_OVERRIDE"):
  mockProgramInOutFilePath=os.environ.get("MOCKPROGRAM_INOUT_FILE_OVERRIDE")
else:
  mockProgramInOutFilePath='.mockprogram_inout.txt'

if not os.path.exists(mockProgramInOutFilePath):
  print("Error: "+mockProgramInOutFilePath+" is missing!")
  sys.exit(1)

mockprogramInout = open(mockProgramInOutFilePath, 'r').read()
mockprogramInoutArray = mockprogramInout.splitlines()
if len(mockprogramInoutArray) and mockprogramInoutArray[-1] == "":
  mockprogramInoutArray = mockprogramInoutArray[:-1]

if len(mockprogramInoutArray) < 3:
  print("Error: "+mockProgramInOutFilePath+" has less than three lines:\n"
        "-------------\n" + mockprogramInout + "-------------")
  sys.exit(2)

# Assert input
expectedInputLine = mockprogramInoutArray[0]
if expectedInputLine.find("MOCK_PROGRAM_INPUT:") != 0:
  print("Error, first line = '" + expectedInputLine + "', does not match "
        "^MOCK_PROGRAM_INPUT:") 
  sys.exit(3)
expectedInput = expectedInputLine.replace("MOCK_PROGRAM_INPUT:", "").strip()
if inputArgs != expectedInput:
  print("Error, input args='" + inputArgs + "' does not match expected='" +
        expectedInput + "'")
  sys.exit(4)

# Get return code
returnCodeLine = mockprogramInoutArray[1]
if returnCodeLine.find("MOCK_PROGRAM_RETURN:") != 0:
  print("Error, second line = '" + returnCodeLine + "', does not match "
        "^MOCK_PROGRAM_RETURN:") 
  sys.exit(5)
returnCode = returnCodeLine.replace("MOCK_PROGRAM_RETURN:", "").strip()

# Get output (can be multi-line)
outputLine = mockprogramInoutArray[2]
if outputLine.find("MOCK_PROGRAM_OUTPUT:") != 0:
  print("Error, third line = '" + outputLine + "', does not match "
        "^MOCK_PROGRAM_OUTPUT:") 
  sys.exit(6)
outputStr = outputLine.replace("MOCK_PROGRAM_OUTPUT: ", "")
numLinesOuput = 1
if len(mockprogramInoutArray) > 3:
  for line in mockprogramInoutArray[3:]:
    if line.find("MOCK_PROGRAM_INPUT:") == 0:
      break
    outputStr = outputStr+"\n"+line
    numLinesOuput = numLinesOuput + 1
print(outputStr)

# Write the remaining lines back into the file
lineLineIndex = 2 + numLinesOuput
if len(mockprogramInoutArray) > lineLineIndex:
  open(mockProgramInOutFilePath, 'w').write(
    ('\n'.join(mockprogramInoutArray[lineLineIndex:]))+"\n" )
else:
  open(mockProgramInOutFilePath, 'w').write("")

# Return exit code
sys.exit(int(returnCode))
