
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


#
# Determine if a given file should be considered a global build file such that
# all of Trilinos should be rebuilt.  In general, any file under the cmake/
# directory with the extension *.cmake is considered a file that requires a
# global rebuild of all Trilinos packaes.  There are a few special files
# that we don't do a global rebuild for:
#
# cmake/TrilinosPackages.cmake: This file gets modified frequently to add new
# packages and rearrange packages.  We don't need to do a global rebuild
# because this list of pakages is validated if we do even a single rebuild.
# If a package line gets removed, the code that reads the Dependencies.cmake
# files will fail and stop.
#
# cmake/TrilinosTPLs.cmake: This file also gets modified frequently.  We
# don't need to enable all Trilinos packages either for the same reason as for
# the TrilinosPackages.cmake file.
#
# EXTRA_REPO/ExtraTrilinosPackages.cmake, EXTRA_REPO/ExtraTrilinosTPLs.cmake
# For external repositories, for the same argument as above, we don't want to
# enable all packages just because these files get modified.
#
# cmake/TPLs/*.cmake: Any FileTPLSOMETHING.cmake file that is not for BLAS,
# LAPACK, or MPI is not needed for Primary Stable code and therefore does not
# need to trigger a global rebulid.
#

def isGlobalBuildFileRequiringGlobalRebuild(modifiedFileFullPath):
  modifiedFileFullPathArray = getFilePathArray(modifiedFileFullPath)
  if len(modifiedFileFullPathArray)==1:
    if modifiedFileFullPathArray[0] == "CMakeLists.txt":
      return True
    if modifiedFileFullPathArray[0] == "Trilinos_version.h":
      return True
  if modifiedFileFullPathArray[0] == 'cmake':
    if modifiedFileFullPathArray[1] == 'TrilinosPackages.cmake':
      return False
    if modifiedFileFullPathArray[1] == 'TrilinosTPLs.cmake':
      return False
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


def getPackageStructFromPath(trilinosDependencies, fullPath):
  packageName = getPackageNameFromPath(trilinosDependencies, fullPath)
  if packageName:
    return trilinosDependencies.getPackageByName(packageName)
  return None


def getPackageNameFromPath(trilinosDependencies, fullPath, prefixPath="packages"):
  return trilinosDependencies.getPackageNameFromPath(fullPath, prefixPath)


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
    if allPackages and isGlobalBuildFileRequiringGlobalRebuild(filePath) and not enabledAllPackages:
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
