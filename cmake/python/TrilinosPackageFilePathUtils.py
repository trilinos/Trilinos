
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


def isGlobalBuildFile(modifiedFileFullPath):
  modifiedFileFullPathArray = getFilePathArray(modifiedFileFullPath)
  if len(modifiedFileFullPathArray)==1:
    if modifiedFileFullPathArray[0] == "CMakeLists.txt":
      return True
    if modifiedFileFullPathArray[0] == "Trilinos_version.h":
      return True
  if modifiedFileFullPathArray[0] == 'cmake':
    if modifiedFileFullPathArray[1] == 'ctest':
      return False
    if modifiedFileFullPathArray[1] == 'DependencyUnitTests':
      return False
    if modifiedFileFullPathArray[1] == 'TPLs':
      if ['FindTPLBLAS.cmake', 'FindTPLLAPACK.cmake', 'FindTPLMPI.cmake'].count(
        modifiedFileFullPathArray[2]) == 1 \
        :
        return True
      return False
    if modifiedFileFullPathArray[-1].rfind(".cmake") != -1:
      return True
  return False


def getFilePathArray(filePathStr):
  return filePathStr.split('/')


def getPackageStructFromPath(trilinosDependencies, fullPath):
  fullPathArray = getFilePathArray(fullPath)
  if fullPathArray[0] == "packages":
    packageDir = fullPathArray[1]
  else:
    packageDir = "../"+fullPathArray[0]
  return trilinosDependencies.getPackageByDir(packageDir)


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


def getPackagesListFromFilePathsList(trilinosDependencies, filePathsList,
  allPackages=False \
  ):
  packagesList = []
  enabledAllPackages = False
  for filePath in filePathsList:
    packageName = getPackageNameFromPath(trilinosDependencies, filePath)
    if findInSequence(packagesList, packageName) == -1 and packageName: 
      packagesList.append(packageName.strip())
    if allPackages and isGlobalBuildFile(filePath) and not enabledAllPackages:
      packagesList.append("ALL_PACKAGES")
      enabledAllPackages = True
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
