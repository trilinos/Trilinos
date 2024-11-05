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
import time


#
# Default logic to determine if a changed file should trigger testing all of
# the packages.
#
# In general, any file at the top level <projectDir>/ or under the
# <projectDir>/cmake/ directory with the extension *.cmake is considered a
# file that requires a global rebuild of all packaes.  However, there are a
# few special files that we don't want to have to do a global rebuild for by
# default.
#
# PackagesList.cmake: This file gets modified frequently to add new packages
# and rearrange packages.  We don't need to do a global rebuild because this
# list of packages is validated if we do even a single rebuild.  If a package
# line gets removed, the code that reads the Dependencies.cmake files will
# fail and stop. However, this has the risk that if only the test test
# category changes (e.g. from 'ST' to 'PT') then this would not trigger a
# rebuild.  But if the package itself was also modified, then that would be
# sufficient for testing.
#
# TPLsList.cmake: This file also gets modified frequently.  We don't need to
# enable all packages either for the same reason as for the PackagesList.cmake
# file.
#
# cmake/ctest/: These drive nightly automated builds and should not require
# testing every package.
#
# cmake/TPLs/*.cmake: Any FileTPLSOMETHING.cmake file is assumed not to be
# needed for Primary Tested code and therefore does not need to trigger a
# global rebuild.
#

class DefaultProjectCiFileChangeLogic:

  def isGlobalBuildFileRequiringGlobalRebuild(self, modifiedFileFullPath):
    modifiedFileFullPathArray = getFilePathArray(modifiedFileFullPath)
    if len(modifiedFileFullPathArray)==1:
      # Files sitting directly under <projectDir>/
      if modifiedFileFullPathArray[0] == "CMakeLists.txt":
        return True
      if modifiedFileFullPathArray[0] == 'PackagesList.cmake':
        return False
      if modifiedFileFullPathArray[0] == 'TPLsList.cmake':
        return False
      if modifiedFileFullPathArray[0].rfind(".cmake") != -1:
        return True
    elif modifiedFileFullPathArray[0] == 'cmake':
      # Files under <projectDir>/cmake/
      if modifiedFileFullPathArray[1]=='ExtraRepositoriesList.cmake':
        return False
      if modifiedFileFullPathArray[1] == 'ctest':
        return False
      if modifiedFileFullPathArray[1] == 'TPLs':
        return False
      if modifiedFileFullPath.rfind("UnitTests/") != -1:
        return False
      if modifiedFileFullPath.rfind(".cmake") != -1:
        return True
    return False


def getProjectCiFileChangeLogic(projectDir):

  if not projectDir:
    return DefaultProjectCiFileChangeLogic()
  else:
    projectCiFileChangeLogicFileBaseDir = os.path.join(projectDir, "cmake")
    projectCiFileChangeLogicFile = os.path.join(
      projectCiFileChangeLogicFileBaseDir, "ProjectCiFileChangeLogic.py")
    if not os.path.isfile(projectCiFileChangeLogicFile):
      return DefaultProjectCiFileChangeLogic()

  # Else, if we get here, then the file
  # <projectDir>/cmake/ProjectCiFileChangeLogic.py exists so lets read it in!

  tribitsDirPath = os.path.abspath(
    os.path.join(
      os.path.dirname(os.path.abspath(__file__)),
      "../..", "tribits"
      )
    )

  old_sys_path = sys.path

  try:
    sys.path = [projectCiFileChangeLogicFileBaseDir] + sys.path
    import ProjectCiFileChangeLogic
    return ProjectCiFileChangeLogic.ProjectCiFileChangeLogic()
  finally:
    sys.path = old_sys_path


def getPackageStructFromPath(projectDependencies, fullPath):
  packageName = getPackageNameFromPath(projectDependencies, fullPath)
  if packageName:
    return projectDependencies.getPackageByName(packageName)
  return None


def getPackageNameFromPath(projectDependencies, fullPath):
  return projectDependencies.getPackageNameFromPath(fullPath)


def extractFilesListMatchingPattern(fileList_in, reMatachingPattern):
  fileList_out = []
  for line in fileList_in:
    reFilePatternMatch = reMatachingPattern.match(line)
    if reFilePatternMatch:
      fileList_out.append(reFilePatternMatch.group(1).strip())
  return fileList_out


def getPackagesListFromFilePathsList(projectDependencies, filePathsList,
  allPackages=False, projectCiFileChangeLogic=DefaultProjectCiFileChangeLogic() \
  ):
  packagesList = []
  enabledAllPackages = False
  for filePath in filePathsList:
    packageName = getPackageNameFromPath(projectDependencies, filePath)
    if findInSequence(packagesList, packageName) == -1 and packageName: 
      packagesList.append(packageName.strip())
    if allPackages \
      and projectCiFileChangeLogic.isGlobalBuildFileRequiringGlobalRebuild(filePath) \
      and not enabledAllPackages \
      :
      packagesList.append("ALL_PACKAGES")
      enabledAllPackages = True
  return packagesList


def getPackageCheckinEmailAddressesListFromFilePathsList(
  projectDependencies, filePathsList \
  ) \
  :
  packageCheckinEmailAddresses = []
  for filePath in filePathsList:
    packageStruct = getPackageStructFromPath(projectDependencies, filePath)
    if packageStruct:
      checkinEmail = packageStruct.emailAddresses.checkin
    else:
      checkinEmail = None
    if findInSequence(packageCheckinEmailAddresses, checkinEmail) == -1:
      packageCheckinEmailAddresses.append(checkinEmail)
  return packageCheckinEmailAddresses
