
#
# ToDo:
#
#  (*) Turn off framework tests by default and turn them in checkin
#      testing ...
#  (*) Turn off generation of HTML/XML files by default and turn them on in
#      checkin testing ...
#

#
# General scripting support
#
# NOTE: Included first to check the version of python!
#

from TrilinosDependencies import getTrilinosDependenciesFromXmlFile
from TrilinosDependencies import defaultTrilinosDepsXmlInFile
from GeneralScriptSupport import *
import time


def getFilePathArray(filePathStr):
  return filePathStr.split('/')


def getPackageStructFromPath(trilinosDependencies, fullPath):
  fullPathArray = getFilePathArray(fullPath)
  if fullPathArray[0] == "packages":
    packageDir = fullPathArray[1]
    return trilinosDependencies.getPackageByDir(packageDir)
  return None


def getPackageNameFromPath(trilinosDependencies, fullPath):
  packageStruct = getPackageStructFromPath(trilinosDependencies, fullPath)
  if packageStruct:
    return packageStruct.packageName
  return u""


def extractFilesListMatchingPattern(fileList_in, reMatachingPattern):
  fileList_out = []
  for line in fileList_in:
    reFilePatternMatch = reMatachingPattern.match(line)
    if reFilePatternMatch:
      fileList_out.append(reFilePatternMatch.group(1).strip())
  return fileList_out


def getPackgesListFromFilePathsList(trilinosDependencies, filePathsList):
  packagesList = []
  for filePath in filePathsList:
    packageName = getPackageNameFromPath(trilinosDependencies, filePath)
    if findInSequence(packagesList, packageName) == -1:
      packagesList.append(packageName)
  return packagesList


def getPackageCheckinEmailAddressesListFromFilePathsList(
  trilinosDependencies, filePathsList \
  ) \
  :
  packageCheckinEmailAddresses = []
  for filePath in filePathsList:
    packageStruct = getPackageStructFromPath(trilinosDependencies, filePath)
    if packageStruct:
      checkinEmail = packageStruct.emailAddresses.checkin
    else:
      # If not a package list then send to all of Trilinos!
      checkinEmail = "trilinos-checkins@software.sandia.gov"
    if findInSequence(packageCheckinEmailAddresses, checkinEmail) == -1:
      packageCheckinEmailAddresses.append(checkinEmail)
  return packageCheckinEmailAddresses


#
# Unit testing code
# 


import unittest


trilinosDependencies = getTrilinosDependenciesFromXmlFile(defaultTrilinosDepsXmlInFile)
#print "\ntrilinosDependencies:\n", trilinosDependencies


updateOutputStr = """
? packages/triutils/doc/html
M cmake/python/checkin-test.py
M cmake/python/dump-cdash-deps-xml-file.py
A packages/nox/src/dummy.C
P packages/stratimikos/dummy.blah
M packages/thyra/src/Thyra_ConfigDefs.hpp
M packages/thyra/CMakeLists.txt
"""

updateOutputList = updateOutputStr.split("\n")


class testTrilinosPackageFilePathUtils(unittest.TestCase):


  def test_getPackageNameFromPath_01(self):
    self.assertEqual(
      getPackageNameFromPath( trilinosDependencies, 'packages/teuchos/CMakeLists.txt' ), 'Teuchos' )


  def test_getPackageNameFromPath_02(self):
    self.assertEqual(
      getPackageNameFromPath( trilinosDependencies, 'packages/thyra/src/blob.cpp' ), 'Thyra' )


  def test_getPackageNameFromPath_noMatch(self):
    self.assertEqual(
      getPackageNameFromPath( trilinosDependencies, 'packages/blob/blob' ), '' )


  def test_extractFilesListMatchingPattern_01(self):

    modifedFilesList = extractFilesListMatchingPattern( updateOutputList,
      re.compile(r"^[MA] (.+)$") )

    modifedFilesList_expected = \
      [
        "cmake/python/checkin-test.py",
        "cmake/python/dump-cdash-deps-xml-file.py",
        "packages/nox/src/dummy.C",
        "packages/thyra/src/Thyra_ConfigDefs.hpp",
        "packages/thyra/CMakeLists.txt",
      ]

    self.assertEqual( modifedFilesList, modifedFilesList_expected )


  def test_getPackgesListFromFilePathsList_01(self):

    filesList = extractFilesListMatchingPattern( updateOutputList,
      re.compile(r"^[AMP?] (.+)$") )
    
    packagesList = getPackgesListFromFilePathsList( trilinosDependencies, filesList )

    packagesList_expected = [u"Triutils", u"", u"NOX", u"Stratimikos", u"Thyra"]

    self.assertEqual( packagesList, packagesList_expected )


  def test_getPackageCheckinEmailAddressesListFromFilePathsList_01(self):

    filesList = extractFilesListMatchingPattern( updateOutputList,
      re.compile(r"^[AMP?] (.+)$") )
    
    packagesList = getPackageCheckinEmailAddressesListFromFilePathsList(
      trilinosDependencies, filesList )

    packagesList_expected = [
      u"triutils-checkins@software.sandia.gov",
      u"trilinos-checkins@software.sandia.gov",
      u"nox-checkins@software.sandia.gov",
      u"stratimikos-checkins@software.sandia.gov",
      u"thyra-checkins@software.sandia.gov"
      ]

    self.assertEqual( packagesList, packagesList_expected )


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(testTrilinosPackageFilePathUtils))
    return suite


if __name__ == '__main__':
  unittest.TextTestRunner(verbosity=2).run(suite())
