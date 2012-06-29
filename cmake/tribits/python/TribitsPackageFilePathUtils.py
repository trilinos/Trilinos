# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
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
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER


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

from TribitsDependencies import getProjectDependenciesFromXmlFile
from GeneralScriptSupport import *
import time


#
# Determine if a given file should be considered a global build file such that
# all of Trilinos should be rebuilt.  In general, any file under the cmake/
# directory with the extension *.cmake is considered a file that requires a
# global rebuild of all Trilinos packaes.  There are a few special files
# that we don't do a global rebuild for:
#
# cmake/ProjectPackages.cmake: This file gets modified frequently to add new
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
    if modifiedFileFullPathArray[0] == 'PackagesList.cmake':
      return False
    if modifiedFileFullPathArray[0] == 'TPLsList.cmake':
      return False
    if modifiedFileFullPathArray[0].rfind(".cmake") != -1:
      return True
  elif modifiedFileFullPathArray[0] == 'cmake':
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


def getPackageStructFromPath(trilinosDependencies, fullPath):
  packageName = getPackageNameFromPath(trilinosDependencies, fullPath)
  if packageName:
    return trilinosDependencies.getPackageByName(packageName)
  return None


def getPackageNameFromPath(trilinosDependencies, fullPath):
  return trilinosDependencies.getPackageNameFromPath(fullPath)


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
