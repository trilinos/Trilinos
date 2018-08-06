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

#
# Usage: mockprogram.py [any arguments]
#
# Mock program that takes input arguments and produces stdout by reading from
# a file .mockprogram_inout.txt in the current directory or the file specified
# by the env var MOCKPROGRAM_INOUT_FILE_OVERRIDE (which can be in any
# diretory).  This script is used to take the place of real commands during a
# test that involves calling commands on the commandline.
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
