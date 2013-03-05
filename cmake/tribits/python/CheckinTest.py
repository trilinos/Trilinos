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
#  (*) Create a TaskStatus class and use it to simplify the logic replacing
#  the simple bools.
#

#
# General scripting support
#
# NOTE: Included first to check the version of python!
#

from GeneralScriptSupport import *

from CheckinTestConstants import *
from TribitsDependencies import getProjectDependenciesFromXmlFile
from TribitsDependencies import getDefaultDepsXmlInFile
from TribitsPackageFilePathUtils import *
import time
import pprint

import re

pp = pprint.PrettyPrinter(indent=4)

# Load some default dependencies for some unit tests
projectDependenciesCache = None
def getDefaultProjectDependenices():
  return projectDependenciesCache


def getGitRepoDir(srcDir, gitRepoName):
  if gitRepoName:
    return srcDir+"/"+gitRepoName
  return srcDir


def getGitRepoFileExt(gitRepoName):
  if gitRepoName:
    return "."+gitRepoName
  return ""


def getCommonConfigFileName():
  return "COMMON.config"


def getProjectDependenciesXmlFileName(projectName):
  return projectName+"PackageDependencies.xml"


def getProjectDependenciesXmlGenerateOutputFileName(projectName):
  return projectName+"PackageDependencies.generate.out"


def getProjectExtraReposPythonOutFile(projectName):
  return projectName+"ExtraRepos.py"


def getTribitsGetExtraReposForCheckinTestOututFile(projectName):
  return projectName+"ExtraRepos.generate.out"


def getBuildSpecificConfigFileName(buildTestCaseName):
  return buildTestCaseName+".config"


def getInitialPullOutputFileName(gitRepoName):
  return "pullInitial"+getGitRepoFileExt(gitRepoName)+".out"


def getInitialExtraPullOutputFileName(gitRepoName):
  return "pullInitialExtra"+getGitRepoFileExt(gitRepoName)+".out"


def getInitialPullSuccessFileName():
  return "pullInitial.success"


def getModifiedFilesOutputFileName(gitRepoName):
  return "modifiedFiles"+getGitRepoFileExt(gitRepoName)+".out"


def getFinalPullOutputFileName(gitRepoName):
  return "pullFinal"+getGitRepoFileExt(gitRepoName)+".out"


def getConfigureOutputFileName():
  return "configure.out"


def getConfigureSuccessFileName():
  return "configure.success"


def getBuildOutputFileName():
  return "make.out"


def getBuildSuccessFileName():
  return "make.success"


def getTestOutputFileName():
  return "ctest.out"


def getTestSuccessFileName():
  return "ctest.success"


def getEmailBodyFileName():
  return "email.out"


def getEmailSuccessFileName():
  return "email.success"


def getFinalCommitBodyFileName(gitRepoName):
  return "commitFinalBody"+getGitRepoFileExt(gitRepoName)+".out"


def getFinalCommitOutputFileName(gitRepoName):
  return "commitFinal"+getGitRepoFileExt(gitRepoName)+".out"


def getCommitStatusEmailBodyFileName():
  return "commitStatusEmailBody.out"


def getPushOutputFileName(gitRepoName):
  return "push"+getGitRepoFileExt(gitRepoName)+".out"


def getExtraCommandOutputFileName():
  return "extraCommand.out"


def getHostname():
  return getCmndOutput("hostname", True)


def getEmailAddressesSpaceString(emailAddressesCommasStr):
  emailAddressesList = emailAddressesCommasStr.split(',')
  return ' '.join(emailAddressesList)


def performAnyBuildTestActions(inOptions):
  if inOptions.doConfigure or inOptions.doBuild \
    or inOptions.doTest or inOptions.doAll or inOptions.localDoAll \
    :
    return True
  return False


def performAnyActions(inOptions):
  if performAnyBuildTestActions(inOptions) or inOptions.doPull:
    return True
  return False


def doGenerateOutputFiles(inOptions):
  return performAnyActions(inOptions)


def doRemoveOutputFiles(inOptions):
  return performAnyActions(inOptions)


def assertEgGitVersionHelper(returnedVersion, expectedVersion):
  if returnedVersion != expectedVersion:
    raise Exception("Error, the installed "+returnedVersion+" does not equal the official "\
      +expectedVersion+"!  To turn this check off, pass in --no-eg-git-version-check.")


def setupAndAssertEgGitVersions(inOptions):

  egWhich = getCmndOutput("which eg", True, False)
  if egWhich == "" or re.match(".+no eg.+", egWhich):
    print "Warning, the eg command is not in your path! ("+egWhich+")"
    setattr(inOptions, "eg", os.path.abspath(inOptions.scriptsDir+"/../common_tools/git/eg"))
    print "Setting to default eg in source tree '"+inOptions.eg+"'!"
  else:
    setattr(inOptions, "eg", "eg")

  egVersionOuput = getCmndOutput(inOptions.eg+" --version", True, False)
  egVersionsList = egVersionOuput.split('\n')

  if inOptions.enableEgGitVersionCheck:
    assertEgGitVersionHelper(egVersionsList[0], "eg version "+g_officialEgVersion)
    assertEgGitVersionHelper(egVersionsList[1], "git version "+g_officialGitVersion)


def assertGitRepoExists(inOptions, gitRepo):
  gitRepoDir = getGitRepoDir(inOptions.srcDir, gitRepo.repoDir)
  if not os.path.os.path.exists(gitRepoDir):
    raise Exception("Error, the specified git repo '"+gitRepo.repoName+"' directory"
      " '"+gitRepoDir+"' does not exist!")


def assertPackageNames(optionName, packagesListStr):
  if not packagesListStr:
    return
  for packageName in packagesListStr.split(','):
    if getDefaultProjectDependenices().packageNameToID(packageName) == -1:
      validPackagesListStr = ""
      for i in range(getDefaultProjectDependenices().numPackages()):
        if validPackagesListStr != "":
          validPackagesListStr += ", "
        validPackagesListStr += getDefaultProjectDependenices().getPackageByID(i).packageName
      raise Exception("Error, invalid package name "+packageName+" in " \
        +optionName+"="+packagesListStr \
        +".  The valid package names include: "+validPackagesListStr)


def assertExtraBuildConfigFiles(extraBuilds):
  if not extraBuilds:
    return
  for extraBuild in extraBuilds.split(','):
    extraBuildConfigFile = extraBuild+".config"
    if not os.path.exists(extraBuildConfigFile):
      raise Exception("Error, the extra build configuration file " \
        +extraBuildConfigFile+" does not exit!")


def getRepoSpaceBranchFromOptionStr(extraPullFrom):
  pullFromRepoList = extraPullFrom.split(':')
  if len(pullFromRepoList) < 2:
    raise Exception("Error, the --extra-pull-from='"+pullFromRepo+"' is not a valid" \
      " <repo>:<branch> specification!")
  repo = ':'.join(pullFromRepoList[0:-1])
  branch = pullFromRepoList[-1]
  return repo + " " + branch


def didSinglePullBringChanges(pullOutFileFullPath):
  pullOutFileStr = readStrFromFile(pullOutFileFullPath)
  #print "\npullOutFileStr:\n" + pullOutFileStr
  alreadyUpToDateIdx = pullOutFileStr.find("Already up-to-date")
  #print "alreadyUpToDateIdx = "+str(alreadyUpToDateIdx)
  return alreadyUpToDateIdx == -1


def executePull(gitRepo, inOptions, baseTestDir, outFile, pullFromRepo=None,
  doRebase=False)\
  :
  cmnd = inOptions.eg+" pull"
  if pullFromRepo:
    repoSpaceBranch = getRepoSpaceBranchFromOptionStr(pullFromRepo)
    print "\nPulling in updates from '"+repoSpaceBranch+"' ...\n"
    cmnd += " " + repoSpaceBranch
  else:
    print "\nPulling in updates from 'origin' ..."
    # NOTE: If you do 'eg pull origin <branch>', then the list of locally
    # modified files will be wrong.  I don't know why this is but if instead
    # you do a raw 'eg pull', then the right list of files shows up.
  if doRebase:
    cmnd += " && "+inOptions.eg+" rebase --against origin/"+inOptions.currentBranch
  outFileFullPath = os.path.join(baseTestDir, outFile)
  (updateRtn, updateTimings) = echoRunSysCmnd( cmnd,
    workingDir=getGitRepoDir(inOptions.srcDir, gitRepo.repoDir),
    outFile=outFileFullPath,
    timeCmnd=True, returnTimeCmnd=True, throwExcept=False
    )
  pullGotChanges = didSinglePullBringChanges(outFileFullPath)
  if pullGotChanges:
    print "\n  Pulled changes from this repo!"
  else:
    print "\n  Did not pull any changes from this repo!"
  return (updateRtn, updateTimings, pullGotChanges)


class Timings:
  def __init__(self):
    self.update = -1.0
    self.configure = -1.0
    self.build = -1.0
    self.test = -1.0
  def deepCopy(self):
    copyTimings = Timings()
    copyTimings.update = self.update
    copyTimings.configure = self.configure
    copyTimings.build = self.build
    copyTimings.test = self.test
    return copyTimings
  def totalTime(self):
    tt = 0.0
    if self.update > 0: tt += self.update
    if self.configure > 0: tt += self.configure
    if self.build > 0: tt += self.build
    if self.test > 0: tt += self.test
    return tt


class GitRepo:
  def __init__(self, repoName, repoDir='', repoType='GIT', repoHasPackages=True):
    self.repoName = repoName
    if repoDir:
      self.repoDir = repoDir
    else:
      self.repoDir = repoName
    self.repoType = repoType
    self.repoHasPackages = repoHasPackages
    self.hasChanges = False
    if (self.repoName and self.repoHasPackages) and (self.repoName != self.repoDir):
      raise Exception("ERROR!  For extra repo '"+repoName+"', if repoHasPackages==True" \
        +" then repoDir must be same as repo name, not '"+repoDir+"'!")
    if self.repoType != 'GIT':
      raise Exception("ERROR!  For extra repo '"+repoName+"', the repo type" \
        +" '"+self.repoType+"' is not supported by the checkin-test.py script, only 'GIT'!")
  def __str__(self):
    return "GitRepo{repoName='"+self.repoName+"'" \
      +", repoDir='"+str(self.repoDir)+"'" \
      +", repoType='"+str(self.repoType)+"'" \
      +", repoHasPackages="+str(self.repoHasPackages) \
      +", hasChanges="+str(self.hasChanges) \
      +"}"
  def __rep__(self):
    return str(self)


def getExtraReposPyFileFromCmakeFile(inOptions, extraReposPythonOutFile, \
  consoleOutputFile = None, verbose=False \
  ):
  if inOptions.extraReposFile == "project":
    extraReposFile = inOptions.srcDir+"/cmake/ExtraRepositoriesList.cmake"
  else:
    extraReposFile = inOptions.extraReposFile
  cmnd = "\""+inOptions.withCmake+"\""+ \
    " -DPROJECT_SOURCE_DIR="+inOptions.srcDir+ \
    " -DTRIBITS_BASE_DIR="+inOptions.tribitsDir+ \
    " -DEXTRA_REPOS_FILE="+extraReposFile+ \
    " -DENABLE_KNOWN_EXTERNAL_REPOS_TYPE="+inOptions.extraReposType+ \
    " -DEXTRA_REPOS="+inOptions.extraRepos+ \
    " -DEXTRA_REPOS_PYTHON_OUT_FILE="+extraReposPythonOutFile
  if inOptions.ignoreMissingExtraRepos:
    cmnd += " -DIGNORE_MISSING_EXTRA_REPOSITORIES=TRUE"
  cmnd += \
    " -P "+inOptions.tribitsDir+"/package_arch/TribitsGetExtraReposForCheckinTest.cmake"
  echoRunSysCmnd(cmnd, throwExcept=True, timeCmnd=True, outFile=consoleOutputFile, \
    verbose=verbose)


def translateExtraReposPyToDictGitRepo(extraReposPyDict):
  repoName = extraReposPyDict['NAME']
  repoDir = extraReposPyDict['DIR']
  repoType = extraReposPyDict['REPOTYPE']
  repoHasPackages = (extraReposPyDict['PACKSTAT'] == 'HASPACKAGES')
  return GitRepo(repoName, repoDir, repoType, repoHasPackages)


class TribitsGitRepos:

  def __init__(self):
    self.reset()
    self.__insertMainRepo()
    self.__initFinalize()

  def initFromCommandlineArguments(self, inOptions, consoleOutputFile=None, verbose=True):

    self.reset()

    self.__insertMainRepo()

    if inOptions.extraRepos!="" and \
      (inOptions.extraReposFile=="" or inOptions.extraReposType=="") \
      :
      # Just use the listed set of extra repos with no checking
      for extraRepoName in inOptions.extraRepos.split(","):
        extraRepo = GitRepo(extraRepoName)
        self.__gitRepoList.append(extraRepo)
    elif inOptions.extraReposFile!="" and inOptions.extraReposType!="":
      # Read in the extra repos from file and assert or ignore missing repos, etc.
      extraReposPythonOutFile = getProjectExtraReposPythonOutFile(inOptions.projectName)
      getExtraReposPyFileFromCmakeFile(inOptions, extraReposPythonOutFile, \
        consoleOutputFile=consoleOutputFile, verbose=verbose)
      extraReposPyTxt = readStrFromFile(extraReposPythonOutFile)
      extraReposPyList = eval(extraReposPyTxt)
      for extraRepoDict in extraReposPyList:
        extraRepo = translateExtraReposPyToDictGitRepo(extraRepoDict)
        self.__gitRepoList.append(extraRepo)

    self.__initFinalize()

  def gitRepoList(self):
    return self.__gitRepoList

  def tribitsExtraRepoNamesList(self):
    return self.__tribitsExtraRepoNamesList

  def numTribitsExtraRepos(self):
    return len(self.__tribitsExtraRepoNamesList)

  def __str__(self):
    strRep = "{\n"
    strRep += "  gitRepoList = " + self.__printReposList(self.__gitRepoList)
    strRep += "  tribitsExtraRepoNamesList = "+str(self.__tribitsExtraRepoNamesList)+"\n"
    strRep += "  }\n"
    return strRep

  def reset(self):
    self.__gitRepoList = []
    self.__tribitsExtraRepoNamesList = []
    return self

  # Private

  def __insertMainRepo(self):
    mainRepo = GitRepo("")
    self.__gitRepoList.append(mainRepo) 

  def __printReposList(self, reposList):
    strRep = "[\n"
    for gitRepo in reposList:
      strRep += ("    " + str(gitRepo) + ",\n")
    strRep += "    ]\n"
    return strRep

  def __initFinalize(self):
    self.__tribitsExtraRepoNamesList = []
    for gitRepo in self.__gitRepoList:
      if gitRepo.repoName and gitRepo.repoHasPackages:
        self.__tribitsExtraRepoNamesList.append(gitRepo.repoName)


def createAndGetProjectDependencies(inOptions, baseTestDir, tribitsGitRepos):

  if tribitsGitRepos.numTribitsExtraRepos() > 0:
    print "\nPulling in packages from extra repos: "+\
       ','.join(tribitsGitRepos.tribitsExtraRepoNamesList())+" ..."
    for gitRepo in tribitsGitRepos.gitRepoList():
      assertGitRepoExists(inOptions, gitRepo)        
    if not inOptions.skipDepsUpdate:
      # There are extra repos so we need to build a new list of Project
      # packages to include the add-on packages.
      cmakeArgumentList = [
        "cmake",
        "-DPROJECT_NAME=%s" % inOptions.projectName,
        cmakeScopedDefine(inOptions.projectName, "TRIBITS_DIR", inOptions.tribitsDir),
        "-DPROJECT_SOURCE_DIR="+inOptions.srcDir,
        cmakeScopedDefine(inOptions.projectName, "OUTPUT_FULL_DEPENDENCY_FILES_IN_DIR",
          baseTestDir),
        cmakeScopedDefine(inOptions.projectName, "EXTRA_REPOSITORIES", "\""+\
          ';'.join(tribitsGitRepos.tribitsExtraRepoNamesList())+"\""),
        "-P %s/cmake/tribits/package_arch/TribitsDumpDepsXmlScript.cmake" % inOptions.srcDir,
        ]
      cmnd = ' '.join(cmakeArgumentList)
      echoRunSysCmnd(cmnd,
        workingDir=baseTestDir,
        outFile=baseTestDir+"/"\
          +getProjectDependenciesXmlGenerateOutputFileName(inOptions.projectName),
        timeCmnd=True)
    else:
      print "\nSkipping update of dependencies XML file on request!"
    projectDepsXmlFile = baseTestDir+"/"\
      +getProjectDependenciesXmlFileName(inOptions.projectName)
  else:
    # No extra repos so you can just use the default list of packages
    projectDepsXmlFile = getDefaultDepsXmlInFile(inOptions.srcDir, inOptions.projectName)

  projectDepsXmlFileOverride = os.environ.get("CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE")
  if projectDepsXmlFileOverride:
    print "\nprojectDepsXmlFileOverride="+projectDepsXmlFileOverride
    projectDepsXmlFile = projectDepsXmlFileOverride

  global projectDependenciesCache
  projectDependenciesCache = getProjectDependenciesFromXmlFile(projectDepsXmlFile)







class BuildTestCase:
  def __init__(self, name, runBuildTestCase, validPackageTypesList,
    isDefaultBuild, skipCaseIfNoChangeFromDefaultEnables,
    extraCMakeOptions, buildIdx \
    ):
    self.name = name
    self.runBuildTestCase = runBuildTestCase
    self.validPackageTypesList = validPackageTypesList
    self.isDefaultBuild = isDefaultBuild
    self.skipCaseIfNoChangeFromDefaultEnables = skipCaseIfNoChangeFromDefaultEnables
    self.extraCMakeOptions = extraCMakeOptions
    self.skippedConfigureDueToNoEnables = False
    self.buildIdx = buildIdx
    self.timings = Timings()


def setBuildTestCaseInList(buildTestCaseList_inout,
  name, runBuildTestCase, validPackageTypesList, isDefaultBuild,
  skipCaseIfNoChangeFromDefaultEnables, extraCMakeOptions \
  ):
  buildTestCaseList_inout.append(
    BuildTestCase(name, runBuildTestCase, validPackageTypesList, isDefaultBuild,
      skipCaseIfNoChangeFromDefaultEnables, extraCMakeOptions,
      len(buildTestCaseList_inout) ) )


def writeDefaultCommonConfigFile():

  commonConfigFileName = getCommonConfigFileName()

  if os.path.exists(commonConfigFileName):

    print "\nThe file "+commonConfigFileName+" already exists!"

  else:

    print "\nCreating a default skeleton file "+commonConfigFileName+" ..."

    commonConfigFileStr = \
      "# Fill in the minimum CMake options that are needed to build and link\n" \
      "# that are common to all builds such as the following:\n" \
      "#\n" \
      "#-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON\n" \
      "#-DBUILD_SHARED:BOOL=ON\n" \
      "#-DTPL_BLAS_LIBRARIES:PATH=/usr/local/libblas.a\n" \
      "#-DTPL_LAPACK_LIBRARIES:PATH=/usr/local/liblapack.a\n" \
      "#\n" \
      "# NOTE: Please do not add any options here that would select what pacakges\n" \
      "# get enabled or disabled.\n"

    writeStrToFile(commonConfigFileName, commonConfigFileStr)


def writeDefaultBuildSpecificConfigFile(buildTestCaseName):

  serialOrMpi = buildTestCaseName.split('_')[0]

  buildSpecificConfigFileName = getBuildSpecificConfigFileName(buildTestCaseName)

  if os.path.exists(buildSpecificConfigFileName):

    print "\nThe file "+buildSpecificConfigFileName+" already exists!"

  else:

    print "\nCreating a default skeleton file "+buildSpecificConfigFileName+" ..."

    buildSpecificConfigFileStr = \
      "# Fill in the minimum CMake options that are needed to build and link\n" \
      "# that are specific to the "+serialOrMpi+" build such as:\n" \
      "#\n" \
      "#-DBUILD_SHARED:BOOL=ON\n"

    if serialOrMpi == "MPI":
      buildSpecificConfigFileStr += \
        "#-DMPI_BASE_DIR:PATH=/usr/lib64/openmpi/1.2.7-gcc\n" \
        "#-DMPI_CXX_COMPILER:PATHNAME=/usr/lib64/openmpi/1.2.7-gcc/mpicxx\n" \
        "#-DMPI_C_COMPILER:PATHNAME=/usr/lib64/openmpi/1.2.7-gcc/mpicc\n" \
        "#-DMPI_Fortran_COMPILER:PATHNAME=/usr/lib64/openmpi/1.2.7-gcc/mpif77\n"
    elif serialOrMpi == "SERIAL":
      buildSpecificConfigFileStr += \
        "#-DCMAKE_CXX_COMPILER:PATHNAME=/usr/local/bin/g++\n" \
        "#-DCMAKE_C_COMPILER:PATHNAME=/usr/local/bin/gcc\n" \
        "#-DCMAKE_Fortran_COMPILER:PATHNAME=/usr/local/bin/gfortran\n"
    else:
      raise Exception("Invalid value for serialOrMpi="+serialOrMpi)
      
    buildSpecificConfigFileStr += \
      "#\n" \
      "# NOTE: Please do not add any options here that would change what pacakges\n" \
      "# or TPLs get enabled or disabled.\n"

    writeStrToFile(buildSpecificConfigFileName, buildSpecificConfigFileStr)


def assertNoIllegalEnables(projectName, fileName, cmakeOption):
  
  reTPlEnable = re.compile(r"-DTPL_ENABLE_.+")
  reProjectEnableOn = re.compile(r"-D%s_ENABLE_[a-zA-Z]+.+=ON" % projectName)

  success = True
  if reTPlEnable.match(cmakeOption):
    print "    ERROR: Illegal TPL enable "+cmakeOption+" in "+fileName+"!"    
    success = False
  elif reProjectEnableOn.match(cmakeOption):
    print "    ERROR: Illegal enable "+cmakeOption+" in "+fileName+"!"    
    success = False
  return success


def readAndAppendCMakeOptions(
    projectName,
    fileName,
    cmakeOptions_inout,
    assertNoIllegalEnablesBool):

  success = True

  if not os.path.exists(fileName):
    return

  print "\nAppending options from "+fileName+":"

  cmakeOptionsFile = open(fileName, 'r')

  for line in cmakeOptionsFile:
    if line[0] != '#':
      cmakeOption = line.strip()
      if cmakeOption == "": continue
      print "  Appending: "+cmakeOption
      if assertNoIllegalEnablesBool:
        if not assertNoIllegalEnables(projectName, fileName, cmakeOption):
          success = False
      cmakeOptions_inout.append(cmakeOption)

  return success


reModifiedFiles = re.compile(r"^[MAD]\t(.+)$")


def getCurrentBranchName(inOptions, baseTestDir):
  branchesStr = getCmndOutput(inOptions.eg+" branch", workingDir=inOptions.srcDir)
  for branchName in branchesStr.split('\n'):
    #print "branchName =", branchName
    if branchName[0] == '*':
      currentBranch = branchName.split(' ')[1]
      #print "currentBranch =", currentBranch
      setattr(inOptions, "currentBranch", currentBranch)
      break


def getCurrentDiffOutput(gitRepo, inOptions, baseTestDir):
  echoRunSysCmnd(
    inOptions.eg+" diff --name-status origin/"+inOptions.currentBranch,
    workingDir=getGitRepoDir(inOptions.srcDir, gitRepo.repoDir),
    outFile=os.path.join(baseTestDir, getModifiedFilesOutputFileName(gitRepo.repoName)),
    timeCmnd=True
    )


def repoHasModifiedFiles(gitRepo, baseTestDir):
  modifiedFilesStr = readStrFromFile(
    baseTestDir+"/"+getModifiedFilesOutputFileName(gitRepo.repoName))
  if modifiedFilesStr:
    return True
  return False


def getCurrentDiffOutputAndLogModified(inOptions, gitRepo, baseTestDir):
  getCurrentDiffOutput(gitRepo, inOptions, baseTestDir)
  gitRepo.hasChanges = repoHasModifiedFiles(gitRepo, baseTestDir)
  if gitRepo.hasChanges:
    print "\n  ==> '"+gitRepo.repoName+"': Has modified files!"
  else:
    print "\n  ==> '"+gitRepo.repoName+"': Does *not* have any modified files!"


def extractPackageEnablesFromChangeStatus(changedFileDiffOutputStr, inOptions_inout,
  gitRepo, enablePackagesList_inout, verbose=True,
  projectDependenciesLocal=None ) \
  :

  if not projectDependenciesLocal:
    projectDependenciesLocal = getDefaultProjectDependenices()

  modifiedFilesList = extractFilesListMatchingPattern(
    changedFileDiffOutputStr.split('\n'), reModifiedFiles )

  for modifiedFileFullPath in modifiedFilesList:

    # Only look for global rebuild files in the master repo (not in extra repos)
    if gitRepo.repoName == '' and isGlobalBuildFileRequiringGlobalRebuild(modifiedFileFullPath):
      if inOptions_inout.enableAllPackages == 'auto':
        if verbose:
          print "\nModifed file: '"+modifiedFileFullPath+"'\n" \
            "  => Enabling all "+inOptions_inout.projectName+" packages!"
        inOptions_inout.enableAllPackages = 'on'

    if gitRepo.repoDir:
      modifiedFileFullPath = gitRepo.repoDir+"/"+modifiedFileFullPath
    #print "\nmodifiedFileFullPath =", modifiedFileFullPath

    packageName = getPackageNameFromPath(projectDependenciesLocal, modifiedFileFullPath)
    if packageName and findInSequence(enablePackagesList_inout, packageName) == -1:
      if verbose:
        print "\nModified file: '"+modifiedFileFullPath+"'\n" \
          "  => Enabling '"+packageName+"'!"
      enablePackagesList_inout.append(packageName)


def createConfigureFile(cmakeOptions, baseCmnd, srcDir, configFileName):

    doConfigStr = ""
  
    doConfigStr += \
      "EXTRA_ARGS=$@\n" \
      "\n" \
      +baseCmnd+ " \\\n"
  
    for opt in cmakeOptions:
      doConfigStr += opt + " \\\n"
    
    doConfigStr += \
      "$EXTRA_ARGS"

    if srcDir:
      doConfigStr += " \\\n"+srcDir
    
    doConfigStr += "\n"
  
    writeStrToFile(configFileName, doConfigStr)
    echoRunSysCmnd('chmod a+x '+configFileName)


def formatMinutesStr(timeInMinutes):
  return ("%.2f" % timeInMinutes) + " min"
  

def getStageStatus(stageName, stageDoBool, stagePassed, stageTiming):
  stageStatusStr = stageName + ": "
  if stageDoBool:
    if stagePassed:
      stageStatusStr += "Passed"
    else:
      stageStatusStr += "FAILED"
    stageStatusStr += " ("+formatMinutesStr(stageTiming)+")"
  else:
    stageStatusStr += "Not Performed"
  stageStatusStr += "\n"
  return stageStatusStr


def getTotalTimeBeginStr(buildTestCaseName):
  return "Total time for "+buildTestCaseName


def getTotalTimeLineStr(buildTestCaseName, timeInMin):
  return getTotalTimeBeginStr(buildTestCaseName)+" = "+formatMinutesStr(timeInMin)


def getTimeInMinFromTotalTimeLine(buildTestCaseName, totalTimeLine):
  if not totalTimeLine:
    return -1.0
  m = re.match(getTotalTimeBeginStr(buildTestCaseName)+r" = (.+) min", totalTimeLine)
  if m and m.groups():
    return float(m.groups()[0])
  else:
    return -1.0


reCtestFailTotal = re.compile(r".+, ([0-9]+) tests failed out of ([0-9]+)")


def analyzeResultsSendEmail(inOptions, buildTestCase,
  enabledPackagesList, cmakeOptions, startingTime, timings ) \
  :

  buildTestCaseName = buildTestCase.name

  print ""
  print "E.1) Determine what passed and failed ..."
  print ""

  success = False

  # Determine if the update passed

  updatePassed = None
  updateOutputExists = False

  if inOptions.doPull:

    if os.path.exists("../"+getInitialPullOutputFileName("")):
      updateOutputExists = True

    if os.path.exists("../"+getInitialPullSuccessFileName()):
      print "\nThe update passed!\n"
      updatePassed = True
    elif updateOutputExists:
      print "\nThe update FAILED!\n"
      updatePassed = False
    else:
      print "\nThe update was never attempted!\n"
      updatePassed = False

  else:

    print "\nThe update step was not performed!\n"

  # Determine if the configured passed

  configurePassed = None
  configureOutputExists = False

  if inOptions.doConfigure:

    if os.path.exists(getConfigureOutputFileName()):
      configureOutputExists = True

    if os.path.exists(getConfigureSuccessFileName()):
      print "\nThe configure passed!\n"
      configurePassed = True
    elif configureOutputExists:
      print "\nThe configure FAILED!\n"
      configurePassed = False
    else:
      print "\nThe configure was never attempted!\n"
      configurePassed = False

  else:

    print "\nThe configure step was not performed!\n"

  # Determine if the build passed

  buildPassed = None
  buildOutputExists = False

  if inOptions.doBuild:

    if os.path.exists(getBuildOutputFileName()):
      buildOutputExists = True

    if os.path.exists(getBuildSuccessFileName()):
      print "\nThe build passed!\n"
      buildPassed = True
    elif buildOutputExists:
      print "\nThe build FAILED!\n"
      buildPassed = False
    else:
      print "\nThe build was never attempted!\n"
      buildPassed = False

  else:

    print "\nThe build step was not performed!\n"

  # Determine if the tests passed

  testsPassed = None
  testOutputExists = False

  if inOptions.doTest:

    if os.path.exists(getTestOutputFileName()):
      testOutputExists = True

    if not testOutputExists:

      print "\nThe tests where never even run!\n"
      testsPassed = False

    else: # testOutputExists

      testResultsLine = getCmndOutput("grep 'tests failed out of' "+getTestOutputFileName(),
        True, False)

      print "testResultsLine = '"+testResultsLine+"'"

      reCtestFailTotalMatch = reCtestFailTotal.match(testResultsLine)

      if reCtestFailTotalMatch:
        numFailedTests = int(reCtestFailTotalMatch.group(1))
        numTotalTests = int(reCtestFailTotalMatch.group(2))
        numPassedTests = numTotalTests - numFailedTests
      else:
        numTotalTests = None
        numPassedTests = None
        testsPassed = False

      if not os.path.exists(getTestSuccessFileName()):
        print "\nThe tests did not run and pass!\n"
        testsPassed = False
      elif numTotalTests == None:
        print "\nCTest was invoked but no tests were run!\n"
        testsPassed = False
      elif numTotalTests == numPassedTests:
        print "\nAll of the tests ran passed!\n"
        testsPassed = True
      else:
        print "\n"+(numTotalTests-numPassedTests)+" tests failed!\n"
        testsPassed = False

  else:

    print "\nRunning the tests was not performed!\n"

  print ""
  print "E.2) Construct the email message ..."
  print ""

  # 2.a) Construct the subject line

  overallPassed = None
  buildCaseStatus = ""
  selectedFinalStatus = False

  if inOptions.doTest and not selectedFinalStatus:
    if testOutputExists:
      if numTotalTests:
        buildCaseStatus += "passed="+str(numPassedTests)+",notpassed="+str(numFailedTests)
      else:
        buildCaseStatus += "no tests run"
      if testsPassed and numTotalTests > 0:
        overallPassed = True
      else:
        overallPassed = False
      selectedFinalStatus = True
    elif not inOptions.doBuild and not buildOutputExists:
      buildCaseStatus += "no active build exists"
      overallPassed = False
      selectedFinalStatus = True

  if inOptions.doBuild and not selectedFinalStatus:
    if buildPassed:
      buildCaseStatus += "build-only passed"
      overallPassed = True
      selectedFinalStatus = True
    elif buildOutputExists:
      buildCaseStatus += "build failed"
      overallPassed = False
      selectedFinalStatus = True

  if inOptions.doConfigure and not selectedFinalStatus:
    if configurePassed:
      buildCaseStatus += "configure-only passed"
      overallPassed = True
      selectedFinalStatus = True
    elif buildTestCase.skippedConfigureDueToNoEnables:
      buildCaseStatus += "skipped configure, build, test due to no enabled packages"
      overallPassed = True
      selectedFinalStatus = True
    elif configureOutputExists:
      buildCaseStatus += "configure failed"
      overallPassed = False
      selectedFinalStatus = True
    else:
      buildCaseStatus += "pre-configure failed"
      overallPassed = False
      selectedFinalStatus = True

  if inOptions.doPull and not selectedFinalStatus:
    if updatePassed:
      buildCaseStatus += "update-only passed"
      overallPassed = True
      selectedFinalStatus = True
    elif updateOutputExists:
      buildCaseStatus += "update FAILED"
      overallPassed = False
      selectedFinalStatus = True

  if not selectedFinalStatus:
    raise Exception("Error, final pass/fail status not found!")

  subjectLine = "%s/%s: %s" % (inOptions.projectName, buildTestCaseName, buildCaseStatus)
  if overallPassed:
    subjectLine = "passed: " + subjectLine
  else:
    subjectLine = "FAILED: " + subjectLine

  print "\nsubjectLine = '"+subjectLine+"'\n"

  success = overallPassed

  # 2.b) Construct the email body

  emailBody = subjectLine + "\n\n"

  emailBody += getCmndOutput("date", True) + "\n\n"

  emailBody += getEnableStatusList(inOptions, enabledPackagesList)
  emailBody += "Hostname: " + getHostname() + "\n"
  emailBody += "Source Dir: " + inOptions.srcDir + "\n"
  emailBody += "Build Dir: " + os.getcwd() + "\n"
  emailBody += "\nCMake Cache Varibles: " + ' '.join(cmakeOptions) + "\n"
  if inOptions.extraCmakeOptions:
    emailBody += "Extra CMake Options: " + inOptions.extraCmakeOptions + "\n"
  if inOptions.makeOptions:
    emailBody += "Make Options: " + inOptions.makeOptions + "\n"
  if inOptions.ctestOptions:
    emailBody += "CTest Options: " + inOptions.ctestOptions + "\n"
  emailBody += "\n"
  emailBody += getStageStatus("Update", inOptions.doPull, updatePassed, timings.update)
  emailBody += getStageStatus("Configure", inOptions.doConfigure, configurePassed, timings.configure)
  emailBody += getStageStatus("Build", inOptions.doBuild, buildPassed, timings.build)
  emailBody += getStageStatus("Test", inOptions.doTest, testsPassed, timings.test)
  emailBody += "\n"

  if inOptions.doTest and testOutputExists and numTotalTests:

    fullCTestOutput = readStrFromFile(getTestOutputFileName())
    if inOptions.showAllTests:
      emailBody += fullCTestOutput
    else:
      emailBody += extractLinesAfterRegex(fullCTestOutput, r".*\% tests passed.*")

  else:

    emailBody += "\n***\n*** WARNING: There are no test results!\n***\n\n"

  endingTime = time.time()
  totalTime = (endingTime - startingTime) / 60.0

  emailBody += "\n"+getTotalTimeLineStr(buildTestCaseName, totalTime)+"\n"

  #print "emailBody:\n\n\n\n", emailBody, "\n\n\n\n"

  writeStrToFile(getEmailBodyFileName(), emailBody)

  if overallPassed:
    echoRunSysCmnd("touch "+getEmailSuccessFileName())

  print ""
  print "E.3) Send the email message ..."
  print ""

  if inOptions.sendEmailTo and buildTestCase.skippedConfigureDueToNoEnables \
     and inOptions.abortGracefullyIfNoEnables \
     :

    print buildTestCaseName + ": Skipping sending build/test case email because" \
      +" there were no enables and --abort-gracefully-if-no-enables was set!"

  elif inOptions.sendEmailTo and inOptions.sendEmailOnlyOnFailure and success:

    print buildTestCaseName + ": Skipping sending build/test case email because" \
     + " it passed and --send-email-only-on-failure was set!"
  
  elif inOptions.sendEmailTo and buildTestCase.skippedConfigureDueToNoEnables \
    and not inOptions.skipCaseSendEmail \
    :

    print "\nSkipping sending final status email for "+buildTestCase.name \
      +" because it had no packages enabled and --skip-case-no-email was set!"

  elif inOptions.sendEmailTo:

    emailAddresses = getEmailAddressesSpaceString(inOptions.sendEmailTo)
    echoRunSysCmnd("mailx -s \""+subjectLine+"\" "+emailAddresses+" < "+getEmailBodyFileName())

  else:

    print "Not sending email because no email addresses were given!"

  # 3) Return final result

  return success


def getBuildTestCaseSummary(testCaseName, trimDown = True):
  # Get the email file
  absEmailBodyFileName = testCaseName+"/"+getEmailBodyFileName()
  if os.path.exists(absEmailBodyFileName):
    testCaseEmailStrArray = open(absEmailBodyFileName, 'r').readlines()
  else:
    testCaseEmailStrArray = None
  # Get the first line (which is the summary)
  testSummaryLine = None
  if testCaseEmailStrArray:
    summaryLine = testCaseEmailStrArray[0].strip()
    if trimDown:
      summaryLineArray = summaryLine.split(":")
      testSummaryLine = summaryLineArray[0].strip() + ": " + summaryLineArray[2].strip()
    else:
      testSummaryLine = summaryLine
  else:
    testSummaryLine = \
      "Error, The build/test was never completed!" \
      " (the file '"+absEmailBodyFileName+"' does not exist.)"
  return testSummaryLine


def getTestCaseEmailSummary(testCaseName, testCaseNum):
  # Get the email file
  absEmailBodyFileName = testCaseName+"/"+getEmailBodyFileName()
  if os.path.exists(absEmailBodyFileName):
    testCaseEmailStrArray = open(absEmailBodyFileName, 'r').readlines()
  else:
    testCaseEmailStrArray = None
  # Write the entry
  testCaseHeader = str(testCaseNum)+") "+testCaseName+" Results:"
  summaryEmailSectionStr = \
    "\n"+testCaseHeader+ \
    "\n"+getStrUnderlineStr(len(testCaseHeader))+"\n" \
    "\n"
  if testCaseEmailStrArray:
    for line in testCaseEmailStrArray:
      summaryEmailSectionStr += "  " + line
    summaryEmailSectionStr += "\n"
  else:
    summaryEmailSectionStr += \
      "Error, The build/test was never completed!" \
      " (the file '"+absEmailBodyFileName+"' does not exist.)\n"
  return summaryEmailSectionStr


def getSummaryEmailSectionStr(inOptions, buildTestCaseList):
  summaryEmailSectionStr = ""
  for buildTestCase in buildTestCaseList:
    if buildTestCase.runBuildTestCase and not buildTestCase.skippedConfigureDueToNoEnables:
      summaryEmailSectionStr += \
        getTestCaseEmailSummary(buildTestCase.name, buildTestCase.buildIdx)
  return summaryEmailSectionStr


def cmakeScopedDefine(projectName, name, value):
  """
  Formats a CMake -D<projectName>_<name>=<value> argument.
  """
  return '-D%s_%s=%s' % (projectName, name, value)


def getEnablesLists(inOptions, validPackageTypesList, isDefaultBuild,
   skipCaseIfNoChangeFromDefaultEnables, tribitsGitRepos,
   baseTestDir, verbose \
   ):

  projectName = inOptions.projectName
  cmakePkgOptions = []
  enablePackagesList = []
  gitRepoList = tribitsGitRepos.gitRepoList()
    
  if inOptions.enablePackages:
    if verbose:
      print "\nEnabling only the explicitly specified packages '"+inOptions.enablePackages+"' ..."
    enablePackagesList = inOptions.enablePackages.split(',')
  else:
    for gitRepo in gitRepoList:
      diffOutFileName = baseTestDir+"/"+getModifiedFilesOutputFileName(gitRepo.repoName)
      if verbose:
        print "\nDetermining the set of packages to enable by examining "+diffOutFileName+" ..."
      if os.path.exists(diffOutFileName):
        changedFileDiffOutputStr = open(diffOutFileName, 'r').read()
        #print "\nchangedFileDiffOutputStr:\n", changedFileDiffOutputStr
        extractPackageEnablesFromChangeStatus(changedFileDiffOutputStr, inOptions, gitRepo,
          enablePackagesList, verbose)
      else:
        if verbose:
          print "\nThe file "+diffOutFileName+" does not exist!\n"

  if verbose:
    print "\nFull package enable list: [" + ','.join(enablePackagesList) + "]"

  if inOptions.disablePackages:
    if verbose:
      print "\nRemoving package enables: [" + inOptions.disablePackages + "]"
    for disablePackage in inOptions.disablePackages.split(","):
      packageIdx = findInSequence(enablePackagesList, disablePackage)
      if packageIdx >= 0:
        del enablePackagesList[packageIdx]

  if verbose:
    print "\nFiltering the set of enabled packages according to allowed package types ..."
  origEnablePackagesList = enablePackagesList[:]
  enablePackagesList = getDefaultProjectDependenices().filterPackageNameList(
    enablePackagesList, validPackageTypesList, verbose)

  if verbose:
    print "\nFinal package enable list: [" + ','.join(enablePackagesList) + "]"

  if tribitsGitRepos.numTribitsExtraRepos() > 0:
    cmakePkgOptions.append(cmakeScopedDefine(
      projectName, "EXTRA_REPOSITORIES:STRING",
      ','.join(tribitsGitRepos.tribitsExtraRepoNamesList())))

  for pkg in enablePackagesList:
    cmakePkgOptions.append(cmakeScopedDefine(projectName, "ENABLE_"+pkg+":BOOL", "ON"))

  cmakePkgOptions.append(cmakeScopedDefine(projectName, "ENABLE_ALL_OPTIONAL_PACKAGES:BOOL", "ON"))

  if inOptions.enableAllPackages == 'on':
    if verbose:
      print "\nEnabling all packages on request!"
    cmakePkgOptions.append(cmakeScopedDefine(projectName, "ENABLE_ALL_PACKAGES:BOOL", "ON"))

  if inOptions.enableFwdPackages:
    if verbose:
      print "\nEnabling forward packages on request!"
    cmakePkgOptions.append(cmakeScopedDefine(projectName, "ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL", "ON"))
  else:
    cmakePkgOptions.append(cmakeScopedDefine(projectName, "ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL", "OFF"))

  if inOptions.disablePackages:
    if verbose:
      print "\nAdding hard disaled for specified packages '"+inOptions.disablePackages+"' ...\n"
    disablePackagesList = inOptions.disablePackages.split(',')
    for pkg in disablePackagesList:
      cmakePkgOptions.append(cmakeScopedDefine(projectName, "ENABLE_"+pkg+":BOOL", "OFF"))

  if verbose:
    print "\ncmakePkgOptions:", cmakePkgOptions

  return (cmakePkgOptions, enablePackagesList)


def runBuildTestCase(inOptions, tribitsGitRepos, buildTestCase, timings):

  success = True

  startingTime = time.time()

  baseTestDir = os.getcwd()

  buildTestCaseName = buildTestCase.name

  if not performAnyActions(inOptions):
    print "\nNo other actions to perform!\n"
    return success

  print "\nCreating a new build directory if it does not already exist ..."
  createDir(buildTestCaseName)

  absBuildDir = os.path.join(baseTestDir, buildTestCaseName)

  echoChDir(absBuildDir)

  try:

    print ""
    print "A) Get the CMake configure options ("+buildTestCaseName+") ..."
    print ""

    preConfigurePassed = True
    projectName = inOptions.projectName

    # A.1) Set the base options
  
    cmakeBaseOptions = []
    if inOptions.extraCmakeOptions:
      cmakeBaseOptions.extend(commandLineOptionsToList(inOptions.extraCmakeOptions))
  
    cmakeBaseOptions.append(cmakeScopedDefine(projectName, "ENABLE_TESTS:BOOL", "ON"))
  
    cmakeBaseOptions.append(cmakeScopedDefine(projectName, "TEST_CATEGORIES:STRING", "BASIC"))

    cmakeBaseOptions.append(cmakeScopedDefine(projectName, "ALLOW_NO_PACKAGES:BOOL", "OFF"))

    if inOptions.ctestTimeOut:
      cmakeBaseOptions.append(("-DDART_TESTING_TIMEOUT:STRING="+str(inOptions.ctestTimeOut)))
  
    cmakeBaseOptions.extend(buildTestCase.extraCMakeOptions)

    result = readAndAppendCMakeOptions(
      inOptions.projectName,
      os.path.join("..", getCommonConfigFileName()),
      cmakeBaseOptions,
      True)
    if not result: preConfigurePassed = False

    reuslt = readAndAppendCMakeOptions(
      inOptions.projectName,
      os.path.join("..", getBuildSpecificConfigFileName(buildTestCaseName)),
      cmakeBaseOptions,
      buildTestCase.isDefaultBuild)
    if not result: preConfigurePassed = False

    print "\ncmakeBaseOptions:", cmakeBaseOptions

    # A.2) Set the package enable options

    cmakePkgOptions = []
    enablePackagesList = []

    if preConfigurePassed:
      (cmakePkgOptions, enablePackagesList) = \
        getEnablesLists(inOptions, buildTestCase.validPackageTypesList, 
          buildTestCase.isDefaultBuild,
          buildTestCase.skipCaseIfNoChangeFromDefaultEnables, tribitsGitRepos,
          baseTestDir, True)
  
    # A.3) Set the combined options

    cmakeOptions = []

    if preConfigurePassed:

      cmakeOptions = cmakeBaseOptions + cmakePkgOptions
    
      print "\ncmakeOptions =", cmakeOptions
    
      print "\nCreating base configure file do-configure.base ..."
      createConfigureFile(cmakeBaseOptions, "cmake", inOptions.srcDir,
        "do-configure.base")
    
      print "\nCreating package-enabled configure file do-configure ..."
      createConfigureFile(cmakePkgOptions, "./do-configure.base", None, "do-configure")
  
    print ""
    print "B) Do the configuration with CMake ("+buildTestCaseName+") ..."
    print ""

    configurePassed = False

    if inOptions.doConfigure and not preConfigurePassed:

      print "\nSKIPPED: "+buildTestCaseName+" configure skipped because pre-configure failed (see above)!\n"

    elif not (enablePackagesList or inOptions.enableAllPackages == 'on'):

      print "\nSKIPPED: "+buildTestCaseName+" configure skipped because no packages are enabled!\n"
      buildTestCase.skippedConfigureDueToNoEnables = True
  
    elif inOptions.doConfigure:
  
      removeIfExists("CMakeCache.txt")
      removeDirIfExists("CMakeFiles")

      cmnd = "./do-configure"

      (configureRtn, timings.configure) = echoRunSysCmnd(cmnd,
        outFile=getConfigureOutputFileName(),
        timeCmnd=True, returnTimeCmnd=True, throwExcept=False
        )

      if configureRtn == 0:
        print "\nConfigure passed!\n"
        echoRunSysCmnd("touch "+getConfigureSuccessFileName())
        configurePassed = True
      else:
        print "\nConfigure failed returning "+str(configureRtn)+"!\n"
        raise Exception("Configure failed!")
  
    else:
  
      print "\nSkipping configure on request!\n"
      if os.path.exists(getConfigureSuccessFileName()):
        print "\nA current successful configure exists!\n"
        configurePassed = True
      else:
        print "\nFAILED: A current successful configure does *not* exist!\n"
  
    print ""
    print "C) Do the build ("+buildTestCaseName+") ..."
    print ""

    buildPassed = False
  
    if inOptions.doBuild and configurePassed:
  
      cmnd = "make"
      if inOptions.makeOptions:
        cmnd += " " + inOptions.makeOptions
  
      (buildRtn, timings.build) = echoRunSysCmnd(cmnd,
        outFile=getBuildOutputFileName(),
        timeCmnd=True, returnTimeCmnd=True, throwExcept=False
        )

      if buildRtn == 0:
        print "\nBuild passed!\n"
        echoRunSysCmnd("touch "+getBuildSuccessFileName())
        buildPassed = True
      else:
        print "\nBuild failed returning "+str(buildRtn)+"!\n"
        raise Exception("Build failed!")
  
    elif inOptions.doBuild and not configurePassed:

      print "\nSKIPPED: "+buildTestCaseName+" build skipped because configure did not pass!\n"
      
    else:

      print "\nSkipping the build on request!\n"
      if os.path.exists(getBuildSuccessFileName()):
        print "\nA current successful build exists!\n"
        buildPassed = True
      else:
        print "\nFAILED: A current successful build does *not* exist!\n"
  
    print ""
    print "D) Run the tests ("+buildTestCaseName+") ..."
    print ""

    testPassed = False

    if inOptions.doTest and buildPassed:
  
      cmnd = "ctest"
      if inOptions.ctestOptions:
        cmnd += " " + inOptions.ctestOptions
  
      (testRtn, timings.test) = echoRunSysCmnd(cmnd,
        outFile=getTestOutputFileName(),
        timeCmnd=True, returnTimeCmnd=True, throwExcept=False
        )
  
      if testRtn == 0:
        print "\nNo tests failed!\n"
        echoRunSysCmnd("touch "+getTestSuccessFileName())
      else:
        errStr = "FAILED: ctest failed returning "+str(testRtn)+"!"
        print "\n"+errStr+"\n"
        raise Exception(errStr)

    elif inOptions.doTest and buildTestCase.skippedConfigureDueToNoEnables:

      print "\nSKIPPED: "+buildTestCaseName+" tests skipped because no packages are enabled!"
      echoRunSysCmnd("touch "+getTestSuccessFileName())
      # NOTE: We have to create this test success file because the presents of
      # this file is used to determine in the build/test case is successful
      # and therefore is okay to push.  This is needed when the script is run
      # a second time to determine if a build/test is successful and therefore
      # allow a push.
  
    else:
  
      print "\nSkipping the tests on request!\n"

  except Exception, e:

    success = False
    printStackTrace()

  print ""
  print "E) Analyze the overall results and send email notification ("+buildTestCaseName+") ..."
  print ""

  if performAnyActions(inOptions):

    result = analyzeResultsSendEmail(inOptions, buildTestCase,
      enablePackagesList, cmakeOptions, startingTime, timings)
    if not result: success = False

  else:

    print "No actions performed, nothing to analyze!"

  return success


def cleanBuildTestCaseOutputFiles(runBuildTestCaseBool, inOptions, baseTestDir, buildTestCaseName):

  if runBuildTestCaseBool and not os.path.exists(buildTestCaseName):

    print "\nSkipping cleaning build/test files for "+buildTestCaseName+" because dir does not exist!\n"

  elif runBuildTestCaseBool and os.path.exists(buildTestCaseName):

    if inOptions.wipeClean:

      print "\nRemoving the existing build directory "+buildTestCaseName+" (--wipe-clean) ..."
      removeDirIfExists(buildTestCaseName)

    elif doRemoveOutputFiles(inOptions):

      echoChDir(buildTestCaseName)
      if inOptions.doConfigure or inOptions.doPull:
        removeIfExists(getConfigureOutputFileName())
        removeIfExists(getConfigureSuccessFileName())
      if inOptions.doBuild or inOptions.doConfigure or inOptions.doPull:
        removeIfExists(getBuildOutputFileName())
        removeIfExists(getBuildSuccessFileName())
      if inOptions.doTest or inOptions.doBuild or inOptions.doConfigure or inOptions.doPull:
        removeIfExists(getTestOutputFileName())
        removeIfExists(getTestSuccessFileName())
      removeIfExists(getEmailBodyFileName())
      removeIfExists(getEmailSuccessFileName())
      echoChDir("..")


def runBuildTestCaseDriver(inOptions, tribitsGitRepos, baseTestDir, buildTestCase, timings):

  success = True

  buildTestCaseName = buildTestCase.name

  print "\n***"
  print "*** Doing build and test of "+buildTestCaseName+" ..."
  print "***\n"
  
  if buildTestCase.runBuildTestCase:

    try:
      echoChDir(baseTestDir)
      writeDefaultBuildSpecificConfigFile(buildTestCaseName)
      result = runBuildTestCase(inOptions, tribitsGitRepos, buildTestCase, timings)
      if not result: success = False
    except Exception, e:
      success = False
      printStackTrace()

  else:

    print "\nSkipping "+buildTestCaseName+" build/test on request!\n"

  return success


def checkBuildTestCaseStatus(buildTestCase, inOptions):

  runBuildTestCaseBool = buildTestCase.runBuildTestCase
  buildTestCaseName = buildTestCase.name
  skippedConfigureDueToNoEnables = buildTestCase.skippedConfigureDueToNoEnables

  statusMsg = None
  timeInMin = -1.0

  if not runBuildTestCaseBool:
    buildTestCaseActionsPass = True
    buildTestCaseOkayToCommit = True
    statusMsg = \
      "Test case "+buildTestCaseName+" was not run! => Does not affect push readiness!"
    return (buildTestCaseActionsPass, buildTestCaseOkayToCommit, statusMsg, timeInMin)

  if skippedConfigureDueToNoEnables:
    buildTestCaseActionsPass = True
    buildTestCaseOkayToCommit = True
    statusMsg = \
      "Skipped configure, build, test due to no enabled packages! => Does not affect push readiness!"
    return (buildTestCaseActionsPass, buildTestCaseOkayToCommit, statusMsg, timeInMin)

  if not os.path.exists(buildTestCaseName) and not performAnyBuildTestActions(inOptions):
    buildTestCaseActionsPass = True
    buildTestCaseOkayToCommit = False
    statusMsg = "No configure, build, or test for "+buildTestCaseName+" was requested!"
    return (buildTestCaseActionsPass, buildTestCaseOkayToCommit, statusMsg, timeInMin)

  if not os.path.exists(buildTestCaseName):
    buildTestCaseActionsPass = False
    buildTestCaseOkayToCommit = False
    statusMsg = "The directory "+buildTestCaseName+" does not exist!"

  emailSuccessFileName = buildTestCaseName+"/"+getEmailSuccessFileName()
  if os.path.exists(emailSuccessFileName):
    buildTestCaseActionsPass = True
  else:
    buildTestCaseActionsPass = False

  testSuccessFileName = buildTestCaseName+"/"+getTestSuccessFileName()
  if os.path.exists(testSuccessFileName):
    buildTestCaseOkayToCommit = True
  else:
    buildTestCaseOkayToCommit = False

  if not statusMsg:
    statusMsg = getBuildTestCaseSummary(buildTestCaseName)

  emailBodyFileName = buildTestCaseName+"/"+getEmailBodyFileName()
  if os.path.exists(emailBodyFileName):
    timeInMinLine = getCmndOutput("grep '"+getTotalTimeBeginStr(buildTestCaseName)+"' " + \
      emailBodyFileName, True, False)
    timeInMin = getTimeInMinFromTotalTimeLine(buildTestCaseName, timeInMinLine)

  return (buildTestCaseActionsPass, buildTestCaseOkayToCommit, statusMsg, timeInMin)


def getUserCommitMessageStr(inOptions):

  absCommitMsgHeaderFile = inOptions.commitMsgHeaderFile
  if not os.path.isabs(absCommitMsgHeaderFile):
    absCommitMsgHeaderFile = os.path.join(inOptions.srcDir, absCommitMsgHeaderFile)

  print "\nExtracting commit message subject and header from the file '" \
        +absCommitMsgHeaderFile+"' ...\n"
  
  commitMsgHeaderFileStr = open(absCommitMsgHeaderFile, 'r').read()
  
  commitEmailBodyStr = commitMsgHeaderFileStr
  
  return commitEmailBodyStr


def getAutomatedStatusSummaryHeaderKeyStr():
  return "Build/Test Cases Summary"


def getAutomatedStatusSummaryHeaderStr():
  
  commitEmailBodyStr = "\n" \
    +getAutomatedStatusSummaryHeaderKeyStr()+"\n"
  
  return commitEmailBodyStr


def getEnableStatusList(inOptions, enabledPackagesList):
  enabledStatusStr = ""
  enabledStatusStr += "Enabled Packages: " + ', '.join(enabledPackagesList) + "\n"
  if inOptions.disablePackages:
    enabledStatusStr += "Disabled Packages: " + inOptions.disablePackages + "\n"
  if inOptions.enableAllPackages == "on":
    enabledStatusStr += "Enabled all Packages\n"
  elif inOptions.enableFwdPackages:
    enabledStatusStr += "Enabled all Forward Packages\n"
  return enabledStatusStr


# Extract the original log message from the output from:
#
#   eg cat-file -p HEAD
#
# This function strips off the git-generated header info and strips off the
# trailing build/test summary data.
#
# NOTE: This function assumes that there will be at least one blank line
# between the buid/test summay data block and the original text message.  If
# there is not, this function will throw!
#
def getLastCommitMessageStrFromRawCommitLogStr(rawLogOutput):

  origLogStrList = []
  pastHeader = False
  numBlankLines = 0
  lastNumBlankLines = 0
  foundStatusHeader = False
  for line in rawLogOutput.split('\n'):
    #print "\nline = '"+line+"'\n"
    if pastHeader:
      origLogStrList.append(line)
      if line == "":
        numBlankLines += 1
      elif numBlankLines > 0:
        lastNumBlankLines = numBlankLines
        numBlankLines = 0
      if line == getAutomatedStatusSummaryHeaderKeyStr():
        foundStatusHeader = True
        break
    if line == "":
      pastHeader = True

  if foundStatusHeader:
    #print "\nlastNumBlankLines =", lastNumBlankLines
    #print "origLogStrList[-1] = '" + origLogStrList[-1] + "'"
    #print "origLogStrList[-2] = '" + origLogStrList[-2] + "'"
    if origLogStrList[-2] != "":
      raise Exception("Error, there must be at least one blank line before the" \
        " build/test summary block!  This is a corrupted commit message.  Please" \
        " use 'git commit --amend' and manually remove the 'Build/test Cases Summary' block.")
    origLogStrList = origLogStrList[0:-lastNumBlankLines]
  else:
    lastNumBlankLines = -1 # Flag we did not find status header

  return ('\n'.join(origLogStrList), lastNumBlankLines)


def getLastCommitMessageStr(inOptions, gitRepo):

  # Get the raw output from the last current commit log
  rawLogOutput = getCmndOutput(
    inOptions.eg+" cat-file -p HEAD",
    workingDir=getGitRepoDir(inOptions.srcDir, gitRepo.repoDir)
    )

  return getLastCommitMessageStrFromRawCommitLogStr(rawLogOutput)[0]


def getLocalCommitsSummariesStr(inOptions, gitRepo, appendRepoName):

  # Get the list of local commits other than this one
  rawLocalCommitsStr = getCmndOutput(
    inOptions.eg+" log --oneline "+inOptions.currentBranch+" ^origin/"+inOptions.currentBranch,
    True,
    workingDir=getGitRepoDir(inOptions.srcDir, gitRepo.repoDir)
    )

  if gitRepo.repoName and appendRepoName:
    repoNameModifier = " ("+gitRepo.repoName+")"
  else:
    repoNameModifier = ""

  print \
    "\nLocal commits for this build/test group"+repoNameModifier+":" \
    "\n----------------------------------------"

  if rawLocalCommitsStr == "\n" or rawLocalCommitsStr == "":
    localCommitsExist = False
  else:
    localCommitsExist = True

  if localCommitsExist:
    print rawLocalCommitsStr
  else:
    print "No local commits exit!"

  localCommitsStr = \
    "Local commits for this build/test group"+repoNameModifier+":\n" \
    "----------------------------------------\n"
  if localCommitsExist:
    localCommitsStr += rawLocalCommitsStr
  else:
    localCommitsStr += "No local commits exist!"

  return (localCommitsStr, localCommitsExist)


def getLocalCommitsSHA1ListStr(inOptions, gitRepo):

  # Get the raw output from the last current commit log
  rawLocalCommitsStr = getCmndOutput(
    inOptions.eg+" log --pretty=format:'%h' "+inOptions.currentBranch+"^ ^origin/"+inOptions.currentBranch,
    True,
    workingDir=getGitRepoDir(inOptions.srcDir, gitRepo.repoDir)
    )

  if rawLocalCommitsStr:
    return ("Other local commits for this build/test group: "
      + (", ".join(rawLocalCommitsStr.split("\n")))) + "\n"

  return ""


def getLocalCommitsExist(inOptions, gitRepo):
  if getLocalCommitsSummariesStr(inOptions, gitRepo, False)[1]:
    return True
  return False


def matchProjectName(line):
  """
  Attempts to match and return the value of PROJECT_NAME in a line like
  SET(PROJECT_NAME <name>)
  If no match can be made, None is returned.
  """
  matchRegex = r'\s*[Ss][Ee][Tt]\s*\(\s*PROJECT_NAME\s+([^\)\s]*)\s*\).*'
  match = re.search(matchRegex, line)
  if match:
    return match.group(1)
  else:
    return None


def getProjectName(sourceDirectory):
  """
  Reads the project name from <root>/ProjectName.cmake
  """
  projectNameFile = os.path.join(sourceDirectory, 'ProjectName.cmake')
  if not os.path.exists(projectNameFile):
    raise Exception(
      "%s is required to exist for a valid Tribits project." % projectNameFile)
  content = open(projectNameFile, "r")
  line = content.readline()
  while line:
    name = matchProjectName(line)
    if name:
      return name
    line = content.readline()
  raise Exception(
    'The file %s does not set the PROJECT_NAME variable. ' +
    'This is required of any Tribits project.')

  
def checkinTest(baseDir, inOptions, configuration={}):
  """
  Main function for checkin testing.
  """

  if not inOptions.projectName:
    inOptions.projectName = getProjectName(inOptions.srcDir)
  
  print "\n**********************************************"
  print "*** Performing checkin testing of %s ***" % inOptions.projectName
  print "**********************************************"

  setattr(inOptions, "tribitsDir", baseDir)

  scriptsDir = os.path.join(baseDir, 'python')
  setattr(inOptions, "scriptsDir", scriptsDir)

  print "\nscriptsDir =", scriptsDir

  print "\nsrcDir =", inOptions.srcDir

  baseTestDir = os.getcwd()
  print "\nbaseTestDir =", baseTestDir

  if inOptions.withoutDefaultBuilds:
    inOptions.defaultBuilds = ''
    
  if inOptions.doAll:
    inOptions.doPull = True
    inOptions.doConfigure = True
    inOptions.doBuild = True
    inOptions.doTest = True

  if inOptions.localDoAll:
    inOptions.allowNoPull = True
    inOptions.doConfigure = True
    inOptions.doBuild = True
    inOptions.doTest = True

  setupAndAssertEgGitVersions(inOptions)

  if inOptions.overallNumProcs:
    inOptions.makeOptions = "-j"+inOptions.overallNumProcs+" "+inOptions.makeOptions
    inOptions.ctestOptions = "-j"+inOptions.overallNumProcs+" "+inOptions.ctestOptions

  assertExtraBuildConfigFiles(inOptions.extraBuilds)
  assertExtraBuildConfigFiles(inOptions.ssExtraBuilds)

  if not inOptions.skipDepsUpdate:
    removeIfExists(getProjectDependenciesXmlFileName(inOptions.projectName))
    removeIfExists(getProjectDependenciesXmlGenerateOutputFileName(inOptions.projectName))
    removeIfExists(getProjectExtraReposPythonOutFile(inOptions.projectName))

  # Set up list of repositories and process dependenices

  tribitsGitRepos = TribitsGitRepos()
  tribitsGitRepos.initFromCommandlineArguments(inOptions)
  #print "\ntribitsGitRepos =", tribitsGitRepos

  createAndGetProjectDependencies(inOptions, baseTestDir, tribitsGitRepos)

  assertPackageNames("--enable-packages", inOptions.enablePackages)
  assertPackageNames("--disable-packages", inOptions.disablePackages)

  success = True

  timings = Timings()

  subjectLine = None
  
  # Set up build/test cases array

  buildTestCaseList = []

  cmakeConfig = configuration.get('cmake', {})

  commonConfigOptions = cmakeConfig.get('common', [])

  defaultBuilds = cmakeConfig.get('default-builds', [])
  requestedDefaultBuilds = inOptions.defaultBuilds
  for buildname, buildopts in defaultBuilds:
    setBuildTestCaseInList(
      buildTestCaseList,
      buildname,
      buildname in requestedDefaultBuilds,
      ["PS"],
      True,
      False,
      commonConfigOptions + buildopts)
      
  if inOptions.ssExtraBuilds:
    for ssExtraBuild in inOptions.ssExtraBuilds.split(','):
      setBuildTestCaseInList(buildTestCaseList, ssExtraBuild, True,
        ["PS", "SS"],  False, True, [])

  allValidPackageTypesList = ["PS", "SS", "EX"]

  if inOptions.extraBuilds:
    for extraBuild in inOptions.extraBuilds.split(','):
      setBuildTestCaseInList(buildTestCaseList, extraBuild, True,
        allValidPackageTypesList,  False, False, [])
  
  try:

    print "\n***"
    print "*** 0) Get the current branch name ..."
    print "***"

    getCurrentBranchName(inOptions, baseTestDir)
    print "\nCurrent branch name = " + inOptions.currentBranch

    print "\n***"
    print "*** 1) Clean old output files ..."
    print "***"

    if inOptions.doPull:
      for gitRepo in tribitsGitRepos.gitRepoList():
        removeIfExists(getInitialPullOutputFileName(gitRepo.repoName))
        removeIfExists(getInitialExtraPullOutputFileName(gitRepo.repoName))
      removeIfExists(getInitialPullSuccessFileName())

    for gitRepo in tribitsGitRepos.gitRepoList():
      removeIfExists(getFinalCommitBodyFileName(gitRepo.repoName))
      removeIfExists(getFinalCommitOutputFileName(gitRepo.repoName))
    removeIfExists(getCommitStatusEmailBodyFileName())

    for gitRepo in tribitsGitRepos.gitRepoList():
      removeIfExists(getModifiedFilesOutputFileName(gitRepo.repoName))
      removeIfExists(getFinalPullOutputFileName(gitRepo.repoName))
      removeIfExists(getPushOutputFileName(gitRepo.repoName))

    if inOptions.executeOnReadyToPush:
      removeIfExists(getExtraCommandOutputFileName())

    for buildTestCase in buildTestCaseList:
      cleanBuildTestCaseOutputFiles(
        buildTestCase.runBuildTestCase, inOptions, baseTestDir, buildTestCase.name)

    print "\n***"
    print "*** 2) Commit changes before pulling updates to merge in (NO LONGER SUPPORTED)"
    print "***"

    print "\n***"
    print "*** 3) Update the %s sources ..." % inOptions.projectName
    print "***"

    repoIsClean = True
    pullPassed = True

    doingAtLeastOnePull = inOptions.doPull

    pulledSomeChanges = False

    if not doingAtLeastOnePull:

      print "\nSkipping all updates on request!\n"

    if doingAtLeastOnePull and pullPassed:

      #
      print "\n3.a) Check that there are no uncommited and no new unknown files before doing the pull(s) ...\n"
      #

      repoIdx = 0
      for gitRepo in tribitsGitRepos.gitRepoList():

        print "\n3.a."+str(repoIdx)+") Git Repo: '"+gitRepo.repoName+"'"

        egStatusOutput = getCmndOutput(inOptions.eg+" status", True, throwOnError=False,
          workingDir=getGitRepoDir(inOptions.srcDir, gitRepo.repoDir))
  
        print \
          "\nOutput from 'eg status':\n" + \
          "\n--------------------------------------------------------------\n" + \
          egStatusOutput + \
          "\n--------------------------------------------------------------\n"

        # See if the repo is clean
  
        if isSubstrInMultiLineString(egStatusOutput, "Changed but not updated"):
          print "\nERROR: There are changed unstaged uncommitted files => cannot continue!"
          repoIsClean = False
  
        if isSubstrInMultiLineString(egStatusOutput, "Changes ready to be committed"):
          print "\nERROR: There are changed staged uncommitted files => cannot continue!"
          repoIsClean = False
  
        if isSubstrInMultiLineString(egStatusOutput, "Newly created unknown files"):
          print "\nERROR: There are newly created uncommitted files => Cannot continue!"
          repoIsClean = False
  
        if not repoIsClean:
          print \
             "\nExplanation: In order to do a meaningful test to allow a push, all files\n" \
             "in the local repo must be committed.  Otherwise, if there are changed but not\n" \
             "committed files or new unknown files that are used in the build or the test, then\n" \
             "what you are testing is *not* what you will be pushing.  If you have changes that\n" \
             "you don't want to push, then try using 'eg stash' before you run this script to\n" \
             "stash away all of the changes you don't want to push.  That way, what you are testing\n" \
             "will be consistent with what you will be pushing.\n"
          pullPassed = False

        #print "gitRepo =", gitRepo
        repoIdx += 1

    if doingAtLeastOnePull and pullPassed:

      # NOTE: We want to pull first from the global repo and then from the
      # extra repo so the extra repo's revisions will get rebased on top of
      # the others.  This is what you would want and expect for the remote
      # test/push process where multiple pulls may be needed before it works.

      #
      print "\n3.b) Pull updates from the global 'origin' repo ..."
      #
    
      if inOptions.doPull and pullPassed:
        repoIdx = 0
        for gitRepo in tribitsGitRepos.gitRepoList():
          print "\n3.b."+str(repoIdx)+") Git Repo: "+gitRepo.repoName
          echoChDir(baseTestDir)
          (updateRtn, updateTimings, pullGotChanges) = executePull(
            gitRepo,
            inOptions, baseTestDir,
            getInitialPullOutputFileName(gitRepo.repoName))
          if pullGotChanges:
            pulledSomeChanges = True
          timings.update += updateTimings
          if updateRtn != 0:
            print "\nPull failed!\n"
            pullPassed = False
            break
          repoIdx += 1
      else:
        print "\nSkipping initial pull from 'origin'!\n"
  
      #
      print "\n3.c) Pull updates from the extra repository '"+inOptions.extraPullFrom+"' ..."
      #

      timings.update = 0
      
      if inOptions.extraPullFrom and pullPassed:
        repoIdx = 0
        for gitRepo in tribitsGitRepos.gitRepoList():
          print "\n3.c."+str(repoIdx)+") Git Repo: "+gitRepo.repoName
          echoChDir(baseTestDir)
          (updateRtn, updateTimings, pullGotChanges) = executePull(
            gitRepo,
            inOptions, baseTestDir,
            getInitialExtraPullOutputFileName(gitRepo.repoName),
            inOptions.extraPullFrom )
          if pullGotChanges:
            pulledSomeChanges = True
          timings.update += updateTimings
          if updateRtn != 0:
            print "\nPull failed!\n"
            pullPassed = False
            break
          repoIdx += 1
      else:
        print "\nSkipping extra pull from '"+inOptions.extraPullFrom+"'!\n"

    # Given overall status of the pulls and determine if to abort gracefully
    if pulledSomeChanges:
      print "There where at least some changes pulled!"
    else:
      print "No changes were pulled!"
 
    #
    print "\nDetermine overall update pass/fail ...\n"
    #

    echoChDir(baseTestDir)

    # Check for prior successful initial pull
    currentSuccessfullPullExists = os.path.exists(getInitialPullSuccessFileName())

    if inOptions.doPull:
      if pullPassed:
        print "\nUpdate passed!\n"
        echoRunSysCmnd("touch "+getInitialPullSuccessFileName())
      else:
        print "\nUpdate failed!\n"
    elif currentSuccessfullPullExists:
      print "\nA previous update was performed and was successful!"
      pullPassed = True
    elif inOptions.allowNoPull:
      print "\nNot performing update since --skip-update was passed in\n"
      pullPassed = True
    else:
      print "\nNo previous successful update is still current!"
      pullPassed = False

    # Update for current successful pull
    currentSuccessfullPullExists = os.path.exists(getInitialPullSuccessFileName())


    print "\n***"
    print "*** 4) Get the list of all the modified files ..."
    print "***"

    hasChangesToPush = False

    if pullPassed:

      for gitRepo in tribitsGitRepos.gitRepoList():
        getCurrentDiffOutputAndLogModified(inOptions, gitRepo, baseTestDir)
        if gitRepo.hasChanges:
          hasChangesToPush = True

    else:

      print "\nSkipping getting list of modified files because pull failed!\n"


    print "\n***"
    print "*** 5) Running the different build/test cases ..."
    print "***"

    # Determine if we will run the build/test cases or not

    # Set runBuildCases flag and other logic
    abortGracefullyDueToNoUpdates = False
    abortGracefullyDueToNoChangesToPush = False
    if not performAnyBuildTestActions(inOptions):
      print "\nNot performing any build cases because no --configure, --build or --test" \
        " was specified!\n"
      runBuildCases = False
    elif doingAtLeastOnePull:
      if repoIsClean and not pulledSomeChanges and \
        inOptions.abortGracefullyIfNoUpdates \
        :
        print "\nNot perfoming any build cases because pull did not bring any commits" \
          " and --abort-gracefully-if-no-updates!\n"
        abortGracefullyDueToNoUpdates = True
        runBuildCases = False
      elif repoIsClean and not hasChangesToPush and \
        inOptions.abortGracefullyIfNoChangesToPush \
        :
        print "\nNot perfoming any build cases because there are no local changes to push" \
          " and --abort-gracefully-if-no-changes-to-push!\n"
        abortGracefullyDueToNoChangesToPush = True
        runBuildCases = False
      elif pullPassed:
        print "\nThe updated passsed, running the build/test cases ...\n"
        runBuildCases = True
      else:
        print "\nNot running any build/test cases because the update (pull) failed!\n"
        runBuildCases = False
    else:
      if inOptions.allowNoPull:
        print "\nNo pull was attemted but we are running the build/test cases anyway" \
          " because --allow-no-pull was specified ...\n"
        runBuildCases = True
      elif os.path.exists(getInitialPullSuccessFileName()):
        print "\nA previous update (pull) was successful, running build/test cases ...!\n"
        runBuildCases = True
      else:
        print "\nNot running any build/test cases because no update was attempted!\n" \
          "\nHint: Use --allow-no-pull to allow build/test cases to run without" \
          " having to do a pull first!"
        runBuildCases = False

    # Run the build/test cases

    buildTestCasesPassed = True

    if runBuildCases:

      echoChDir(baseTestDir)
  
      writeDefaultCommonConfigFile()

      print "\nSetting up to run the build/test cases:"
      for i in range(len(buildTestCaseList)):
        buildTestCase = buildTestCaseList[i]
        print str(i) + ") " + buildTestCase.name + ":",
        if buildTestCase.runBuildTestCase:
          print "Will attempt to run!"
        else:
          print "Will *not* attempt to run on request!"

      for buildTestCase in buildTestCaseList:
        buildTestCase.timings = timings.deepCopy()
        result = runBuildTestCaseDriver(
          inOptions,
          tribitsGitRepos,
          baseTestDir,
          buildTestCase,
          buildTestCase.timings
          )
        if not result:
          buildTestCasesPassed = False
          success = False


    print "\n***"
    print "*** 6) Determine overall success and push readiness ..."
    print "***"

    okayToCommit = False
    okayToPush = False
    forcedCommitPush = False
    abortedCommitPush = False
    atLeastOneConfigureBuildAttemptPassed = False

    if inOptions.doPushReadinessCheck:

      echoChDir(baseTestDir)
  
      okayToCommit = success
      subjectLine = None

      commitEmailBodyExtra = ""
      shortCommitEmailBodyExtra = ""

      (cmakePkgOptions, enabledPackagesList) = \
        getEnablesLists(inOptions, allValidPackageTypesList, False, False,
        tribitsGitRepos, baseTestDir, False)

      enableStatsListStr = getEnableStatusList(inOptions, enabledPackagesList)
      commitEmailBodyExtra += enableStatsListStr
      shortCommitEmailBodyExtra += enableStatsListStr

      commitEmailBodyExtra += \
        "\nBuild test results:" \
        "\n-------------------\n"

      for i in range(len(buildTestCaseList)):
        buildTestCase = buildTestCaseList[i]
        buildTestCaseName = buildTestCase.name
        (buildTestCaseActionsPass, buildTestCaseOkayToCommit, statusMsg, timeInMin) = \
          checkBuildTestCaseStatus(buildTestCase, inOptions)
        buildTestCaseStatusStr = str(i)+") "+buildTestCaseName+" => "+statusMsg
        if not buildTestCaseOkayToCommit:
          buildTestCaseStatusStr += " => Not ready to push!"
        buildTestCaseStatusStr += " ("+formatMinutesStr(timeInMin)+")\n"
        print buildTestCaseStatusStr
        commitEmailBodyExtra += buildTestCaseStatusStr
        shortCommitEmailBodyExtra += buildTestCaseStatusStr
        #print "buildTestCaseActionsPass =", buildTestCaseActionsPass
        if not buildTestCaseActionsPass:
          success = False
        if not buildTestCaseOkayToCommit:
          okayToCommit = False
        #print "buildTestCaseOkayToCommit =", buildTestCaseOkayToCommit
        if buildTestCase.runBuildTestCase and buildTestCaseOkayToCommit \
           and not buildTestCase.skippedConfigureDueToNoEnables \
           :
           #print "Setting atLeastOneConfigureBuildAttemptPassed=True"
           atLeastOneConfigureBuildAttemptPassed = True

      if not atLeastOneConfigureBuildAttemptPassed:
        print "\nThere were no successfuly attempts to configure/build/test!"
        okayToCommit = False

      if not okayToCommit:
        print "\nAt least one of the actions (update, configure, built, test)" \
          " failed or was not performed correctly!\n"
     
      # Determine if we should do a forced push
      if inOptions.doPushReadinessCheck and not okayToCommit and inOptions.forcePush \
        :
        forcedPushMsg = \
          "\n***" \
          "\n*** WARNING: The acceptance criteria for doing a push has *not*" \
          "\n*** been met, but a push is being forced anyway by --force-push!" \
          "\n***\n"
        print forcedPushMsg
        okayToCommit = True
        forcedCommitPush = True

      # Determine if a push is ready to try or not

      if okayToCommit:
        if currentSuccessfullPullExists:
          print "\nA current successful pull also exists => Ready for final push!\n"
          okayToPush = True
        else:
          commitEmailBodyExtra += \
             "\nA current successful pull does *not* exist => Not ready for final push!\n" \
             "\nExplanation: In order to safely push, the local working directory needs\n" \
             "to be up-to-date with the global repo or a full integration has not been\n" \
             "performed!\n"
          print commitEmailBodyExtra
          okayToPush = False
          abortedCommitPush = True
      else:
        okayToPush = False
  
      if okayToPush:
        print "\n  => A PUSH IS READY TO BE PERFORMED!"
      else:
        print "\n  => A PUSH IS *NOT* READY TO BE PERFORMED!"

    else:

      print "\nSkipping push readiness check on request!"
      okayToCommit = False
      okayToPush = False

  
    print "\n***"
    print "*** 7) Do final push  ..."
    print "***"

    # Attempt the final pull, commit amend, and push

    pullFinalPassed = True
    amendFinalCommitPassed = True
    pushPassed = True
    didPush = False
    allLocalCommitSummariesStr = ""
    
    if not inOptions.doPush:
  
      print "\nNot doing the push but sending an email" \
            " about the commit/push readiness status ..."
  
      if okayToPush:
        subjectLine = "READY TO PUSH"
      else:
        subjectLine = "NOT READY TO PUSH"

    elif not okayToPush:

      print "\nNot performing push due to prior errors\n"
      pushPassed = False

    else: # inOptions.doPush and okayToPush:

      #
      print "\n7.a) Performing a final pull to make sure there are no conflicts for push ...\n"
      #
      
      if not okayToPush:

        print "\nSkippng final pull due to prior errors!\n"
        pullFinalPassed = False

      else:

        print "\nExplanation: In order to push, the local repo needs to be up-to-date\n" \
          "with the global repo or the push will not be allowed.  Therefore, a pull\n" \
          "before the push must be performed if there are updates in the global reop\n" \
          "regardless if --pull was specified or not.  Also, a rebase might be done in\n" \
          "order to get a linear history required by the hooks in the main repository.\n"

        doFinalRebase = inOptions.rebase
        if not doFinalRebase:
          print "Skipping the final rebase on request! (see --no-rebase option)"

        pullFinalPassed = True

        repoIdx = 0
        for gitRepo in tribitsGitRepos.gitRepoList():
  
          print "\n7.a."+str(repoIdx)+") Git Repo: '"+gitRepo.repoName+"'"
  
          (update2Rtn, update2Time, pullGotChanges) = \
            executePull(gitRepo, inOptions, baseTestDir,
              getFinalPullOutputFileName(gitRepo.repoName), None,
              doFinalRebase )
  
          if update2Rtn != 0:
            pullFinalPassed = False
            break

          repoIdx += 1

        if pullFinalPassed:
          print "\nFinal update passed!\n"
        else:
          print "\nFinal update failed!\n"

        if not pullFinalPassed: okayToPush = False

      #
      print "\n7.b) Amending the final commit message by appending test results ...\n"
      #

      if not inOptions.appendTestResults:

        print "\nSkipping appending test results on request (--no-append-test-results)!\n"

      elif not okayToPush:

        print "\nSkippng appending test results due to prior errors!\n"
        amendFinalCommitPassed = False

      else:  # inOptions.appendTestResults and okayToPush
  
        print "\nAttempting to amend the final commmit message ...\n"

        repoIdx = 0
        for gitRepo in tribitsGitRepos.gitRepoList():
  
          print "\n7.b."+str(repoIdx)+") Git Repo: '"+gitRepo.repoName+"'"

          try:

            lastCommitMessageStr = getLastCommitMessageStr(inOptions, gitRepo)
            #print "\nlastCommitMessageStr:\n-------------\n"+lastCommitMessageStr+"-------------\n"
            (localCommitSummariesStr, localCommitsExist) = \
              getLocalCommitsSummariesStr(inOptions, gitRepo, False)
            #print "\nlocalCommitsExist =", localCommitsExist, "\n"
            localCommitSHA1ListStr = getLocalCommitsSHA1ListStr(inOptions, gitRepo)
  
            # Get then final commit message
            finalCommitEmailBodyStr = lastCommitMessageStr
            finalCommitEmailBodyStr += getAutomatedStatusSummaryHeaderStr()
            finalCommitEmailBodyStr += shortCommitEmailBodyExtra
            finalCommitEmailBodyStr += localCommitSHA1ListStr
            if forcedCommitPush:
              finalCommitEmailBodyStr += "WARNING: Forced the push!\n"
            finalCommitEmailBodyFileName = getFinalCommitBodyFileName(gitRepo.repoName)
            writeStrToFile(finalCommitEmailBodyFileName, finalCommitEmailBodyStr)
  
            # Amend the final commit message
            if localCommitsExist:
  
              commitAmendRtn = echoRunSysCmnd(
                inOptions.eg+" commit --amend" \
                " -F "+os.path.join(baseTestDir, finalCommitEmailBodyFileName),
                workingDir=getGitRepoDir(inOptions.srcDir, gitRepo.repoDir),
                outFile=os.path.join(baseTestDir, getFinalCommitOutputFileName(gitRepo.repoName)),
                timeCmnd=True, throwExcept=False
                )
  
              if commitAmendRtn != 0:
                amendFinalCommitPassed = False
                break
  
            else:
  
              print "\nSkipping amending last commit because there are no local commits!\n"
              
          except Exception, e:
            success = False
            amendFinalCommitPassed = False
            printStackTrace()

          repoIdx += 1

        # end for

        if amendFinalCommitPassed:
          print "\nAppending test results to last commit passed!\n"
        else:
          print "\nAppending test results to last commit failed!\n"

      if not amendFinalCommitPassed: okayToPush = False

      # Get the updated SHA1 after the commit has been (or has not been)
      # amended but before the push!  NOTE: We grab the list of commits even
      # if we don't ammend the last commit message
      repoIdx = 0
      for gitRepo in tribitsGitRepos.gitRepoList():
        localCommitSummariesStr = \
          getLocalCommitsSummariesStr(inOptions, gitRepo, True)[0]
        if allLocalCommitSummariesStr:
          allLocalCommitSummariesStr += ("\n\n" + localCommitSummariesStr)
        else:
          allLocalCommitSummariesStr = localCommitSummariesStr
        repoIdx += 1

      #
      print "\n7.c) Pushing the the local commits to the global repo ...\n"
      #

      if not okayToPush:

        print "\nNot performing push due to prior errors!\n"
        pushPassed = False

      else:
  
        print "\nAttempting to do the push ..."

        debugSkipPush = os.environ.get("CHECKIN_TEST_SKIP_PUSH","")
        #print "debugSkipPush =", debugSkipPush
        #debugSkipPush = True

        didAtLeastOnePush = False

        repoIdx = 0
        for gitRepo in tribitsGitRepos.gitRepoList():
  
          print "\n7.c."+str(repoIdx)+") Git Repo: '"+gitRepo.repoName+"'"

          if gitRepo.hasChanges:

            if not debugSkipPush:
              pushRtn = echoRunSysCmnd(
                inOptions.eg+" push origin "+inOptions.currentBranch,
                workingDir=getGitRepoDir(inOptions.srcDir, gitRepo.repoDir),
                outFile=os.path.join(baseTestDir, getPushOutputFileName(gitRepo.repoName)),
                throwExcept=False, timeCmnd=True )
              didAtLeastOnePush = True
            else:
              print "\nSkipping push due to debug override ..."
              pushRtn = 0
    
            if pushRtn != 0:
              pushPassed = False
              break

          else:

            print "\nSkipping push to '"+gitRepo.repoName+"' because there are no changes!"
  
          repoIdx += 1

        # end for
  
        if pushPassed:
          if didAtLeastOnePush:
            print "\nPush passed!\n"
            didPush = True
          else:
            print "\nPush failed because the push was never attempted!"
        else:
          print "\nPush failed!\n"

      if not pushPassed: okayToPush = False

  
    print "\n***"
    print "*** 8) Set up to run execute extra command on ready to push  ..."
    print "***"

    if inOptions.executeOnReadyToPush and not okayToPush:

      print "\nNot executing final command ("+inOptions.executeOnReadyToPush+") since" \
        +" a push is not okay to be performed!\n"

    elif inOptions.executeOnReadyToPush and okayToPush:

      executeCmndStr = "\nExecuting final command ("+inOptions.executeOnReadyToPush+") since" \
        +" a push is okay to be performed!\n"
      commitEmailBodyExtra += executeCmndStr
      print executeCmndStr

    else:

      print "\nNot executing final command since none was given ...\n"

  
    print "\n***"
    print "*** 9) Create and send push (or readiness status) notification email  ..."
    print "***\n"
    
    allConfiguresAbortedDueToNoEnablesGracefullAbort = False

    if inOptions.doPushReadinessCheck:

      #
      print "\n9.a) Getting final status to send out in the summary email ...\n"
      #

      # Determine if all configures were aborted because no package enables
      allConfiguresAbortedDueToNoEnablesGracefullAbort = True
      for buildTestCase in buildTestCaseList:
        if not buildTestCase.skippedConfigureDueToNoEnables:
          allConfiguresAbortedDueToNoEnablesGracefullAbort = False

      if not pullPassed:
        subjectLine = "INITIAL PULL FAILED"
        commitEmailBodyExtra += "\n\nFailed because initial pull failed!" \
          " See '"+getInitialPullOutputFileName("*")+"'\n\n"
        success = False
      elif abortGracefullyDueToNoUpdates:
        subjectLine = "ABORTED DUE TO NO UPDATES"
        commitEmailBodyExtra += "\n\nAborted because no updates and --abort-gracefully-if-no-updates was set!\n\n"
        success = True
      elif abortGracefullyDueToNoChangesToPush:
        subjectLine = "ABORTED DUE TO NO CHANGES TO PUSH"
        commitEmailBodyExtra += "\n\nAborted because no changes to push and --abort-gracefully-if-no-changes-to-push was set!\n\n"
        success = True
      elif allConfiguresAbortedDueToNoEnablesGracefullAbort:
        subjectLine = "ABORTED DUE TO NO ENABLES"
        commitEmailBodyExtra += "\n\nAborted because no enables and --abort-gracefully-if-no-enables was set!\n\n"
        success = True
      elif not pullFinalPassed:
        subjectLine = "FINAL PULL FAILED"
        commitEmailBodyExtra += "\n\nFailed because the final pull failed!" \
          " See '"+getFinalPullOutputFileName("*")+"'\n\n"
        success = False
      elif not amendFinalCommitPassed:
        subjectLine = "AMEND COMMIT FAILED"
        commitEmailBodyExtra += "\n\nFailed because the final test commit amend failed!" \
          " See '"+getFinalCommitOutputFileName("*")+"'\n\n"
        success = False
      elif inOptions.doPush and pushPassed and forcedCommitPush:
        subjectLine = "DID FORCED PUSH"
        commitEmailBodyExtra += forcedPushMsg
        success = True
        commitEmailBodyExtra += forcedPushMsg
      elif not buildTestCasesPassed:
        subjectLine = "FAILED CONFIGURE/BUILD/TEST"
        commitEmailBodyExtra += "\n\nFailed because one of the build/test cases failed!\n"
        success = False
      elif inOptions.doPush:
        if didPush and not forcedCommitPush:
          subjectLine = "DID PUSH"
        elif abortedCommitPush:
          subjectLine = "ABORTED COMMIT/PUSH"
          commitEmailBodyExtra += "\n\nCommit/push was never attempted since commit/push" \
          " criteria failed!\n\n"
          success = False
        else:
          subjectLine = "PUSH FAILED"
          commitEmailBodyExtra += "\n\nFailed because push failed!" \
            " See '"+getPushOutputFileName("*")+"'\n\n"
          success = False
      else:
        if okayToPush:
          subjectLine = "READY TO PUSH"
        else:
          subjectLine = "NOT READY TO PUSH"

      #
      print "\n9.b) Create and send out push (or readiness status) notification email ..."
      #
    
      subjectLine += ": %s: %s" % (inOptions.projectName, getHostname())
    
      emailBodyStr = subjectLine + "\n\n"
      emailBodyStr += getCmndOutput("date", True) + "\n\n"
      emailBodyStr += commitEmailBodyExtra + "\n"
      emailBodyStr += allLocalCommitSummariesStr + "\n"
      emailBodyStr += getSummaryEmailSectionStr(inOptions, buildTestCaseList)
    
      print "\nCommit status email being sent:\n" \
        "--------------------------------\n\n\n\n"+emailBodyStr+"\n\n\n\n"
  
      summaryCommitEmailBodyFileName = getCommitStatusEmailBodyFileName()
    
      writeStrToFile(summaryCommitEmailBodyFileName, emailBodyStr)

      if inOptions.sendEmailTo and abortGracefullyDueToNoUpdates:

        print "\nSkipping sending final email because there were no updates" \
          " and --abort-gracefully-if-no-updates was set!"

      elif inOptions.sendEmailTo and abortGracefullyDueToNoChangesToPush:

        print "\nSkipping sending final email because there are no local changes to push" \
          " and --abort-gracefully-if-no-changes-to-push was set!"

      elif inOptions.sendEmailTo and allConfiguresAbortedDueToNoEnablesGracefullAbort:

        print "\nSkipping sending final email because there were no enables" \
          " and --abort-gracefully-if-no-enables was set!"

      elif inOptions.sendEmailTo and inOptions.sendEmailOnlyOnFailure and success:

        print "\nSkipping sending final email because it passed" \
          " and --send-email-only-on-failure was set!"
 
      elif inOptions.sendEmailTo:
  
        emailAddresses = getEmailAddressesSpaceString(inOptions.sendEmailTo)
        if inOptions.sendEmailToOnPush and didPush:
          emailAddresses += " " + getEmailAddressesSpaceString(inOptions.sendEmailToOnPush)
        echoRunSysCmnd("mailx -s \""+subjectLine+"\" " \
          +emailAddresses+" < "+summaryCommitEmailBodyFileName)
  
      else:
  
        print "\nNot sending push readiness status email because --send-email-to is empty!"

    else:

      print "\nNot performing push or sending out push readiness status on request!"

  
    print "\n***"
    print "*** 10) Run execute extra command on ready to push  ..."
    print "***"

    if inOptions.executeOnReadyToPush and okayToPush:

      print executeCmndStr

      extraCommandRtn = echoRunSysCmnd(
        inOptions.executeOnReadyToPush,
        workingDir=baseTestDir,
        outFile=os.path.join(baseTestDir, getExtraCommandOutputFileName()),
        throwExcept=False, timeCmnd=True )

      if extraCommandRtn == 0:
        print "\nExtra command passed!\n"
      else:
        print "\nExtra command failed!\n"
        success = False

    else:

      print "\nNot executing final command ...\n"

  
    if not performAnyActions(inOptions) and not inOptions.doPush:

      print \
        "\n***\n" \
        "*** WARNING: No actions were performed!\n" \
        "***\n" \
        "*** Hint: Specify --do-all to perform full integration update/build/test\n" \
        "*** or --push to push the commits for a previously run test!\n" \
        "***\n\n"
  
  except Exception, e:

    success = False
    printStackTrace()

  g_sysCmndInterceptor.assertAllCommandsRun()

  # Print the final status at the very end
  if subjectLine:
    print \
      "\n\n" + subjectLine + "\n\n"
  
  return success
