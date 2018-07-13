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
# General scripting support
#
# NOTE: Included first to check the version of python!
#

from TribitsDependencies import getProjectDependenciesFromXmlFile
from GeneralScriptSupport import *


def getPackageNameFromTestName(trilinosDependencies, testName):
  return trilinosDependencies.getPackageNameFromTestName(testName)


def getTestNameFromLastTestsFailedLine(trilinosDependencies, line):
  lineArray = line.split(':')
  assert len(lineArray) == 2, "Error, the line '"+line+"' not formatted correctly!"
  testName = lineArray[1]
  assert testName != "", "Error, test name '"+testName+"' can't be empty!"
  return testName


#
# Given the lines from a LastTestsFail*.log file, return an array of the
# matching parent package names.
#
# This will return the list of matching packages only once per package.
#

def getPackageNamesFromLastTestsFailedLines(trilinosDependencies, \
  lastTestsFailedLines \
  ):
  #print ("\nlastTestsFailedLine:\n"+str(lastTestsFailedLines))
  packageNames = []
  for lastTestsFailedLine in lastTestsFailedLines:
    #print ("\nlastTestsFailedLine = '"+lastTestsFailedLine+"'")
    testName = \
      getTestNameFromLastTestsFailedLine(trilinosDependencies, lastTestsFailedLine)
    #print ("\ntestName = '"+testName+"'")
    packageName = getPackageNameFromTestName(trilinosDependencies, testName)
    #print("\npackageName = '"+packageName+"'")
    if findInSequence(packageNames, packageName) == -1 and packageName:
      #print ("\nAppend '"+packageName+"'")
      packageNames.append(packageName)
  return packageNames

