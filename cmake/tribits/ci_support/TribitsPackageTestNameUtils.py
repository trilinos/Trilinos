# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
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

