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


########################################
# Unit testing code for CheckinTest.py #
########################################


from CheckinTest import *
import unittest

scriptsDir = getScriptBaseDir()

tribitsBaseDir=os.path.abspath(scriptsDir+"/..")
packageArchDir=os.path.abspath(tribitsBaseDir+"/package_arch")
mockProjectBaseDir=os.path.abspath(packageArchDir+"/UnitTests/MockTrilinos")

#####################################################
#
# Testing helper code
#
#####################################################


class MockOptions:
  def __init__(self):
    self.projectName = "Trilinos"
    self.srcDir = mockProjectBaseDir
    self.tribitsDir = tribitsBaseDir 
    self.enableAllPackages = 'auto'
    self.extraReposFile = ""
    self.extraReposType = ""
    self.extraRepos = ""
    self.ignoreMissingExtraRepos = ""
    self.withCmake = "cmake"


def assertGrepFileForRegexStrList(testObject, testName, fileName, regexStrList, verbose):
  assert(os.path.isfile(fileName))
  for regexToFind in regexStrList.strip().split('\n'):
    if regexToFind == "": continue
    foundRegex = getCmndOutput("grep '"+regexToFind+"' "+fileName, True, False)
    if verbose or not foundRegex:
      print "\n"+testName+": In '"+fileName+"' look for regex '"+regexToFind+"' ...", 
      print "'"+foundRegex+"'", 
      if foundRegex: print ": PASSED"
      else: print ": FAILED"
    testObject.assertNotEqual(foundRegex, "")


def assertNotGrepFileForRegexStrList(testObject, testName, fileName, regexStrList, verbose):
  assert(os.path.isfile(fileName))
  for regexToFind in regexStrList.strip().split('\n'):
    if regexToFind == "": continue
    foundRegex = getCmndOutput("grep '"+regexToFind+"' "+fileName, True, False)
    if verbose or foundRegex:
      print "\n"+testName+": In '"+fileName \
        +"' assert not exist regex '"+regexToFind+"' ... '"+foundRegex+"'", 
      if foundRegex: print ": FAILED"
      else: print ": PASSED"
    testObject.assertEqual(foundRegex, "")


#############################################################################
#
# Test formatMinutesStr
#
#############################################################################


class test_formatMinutesStr(unittest.TestCase):

  def test_00(self):
    self.assertEqual(formatMinutesStr(0.000000), "0.00 min")

  def test_01(self):
    self.assertEqual(formatMinutesStr(1245.244678), "1245.24 min")

  def test_02(self):
    self.assertEqual(formatMinutesStr(1245.245678), "1245.25 min")

  def test_03(self):
    self.assertEqual(formatMinutesStr(1.245678), "1.25 min")

  def test_04(self):
    self.assertEqual(formatMinutesStr(0.202), "0.20 min")

  def test_05(self):
    self.assertEqual(formatMinutesStr(0.204), "0.20 min")

  def test_06(self):
    self.assertEqual(formatMinutesStr(0.2053333), "0.21 min")

  def test_07(self):
    self.assertEqual(formatMinutesStr(0.2943333), "0.29 min")

  def test_08(self):
    self.assertEqual(formatMinutesStr(0.2993333), "0.30 min")

  def test_09(self):
    self.assertEqual(formatMinutesStr(45.2993333), "45.30 min")

  def test_10(self):
    self.assertEqual(formatMinutesStr(45.2493333), "45.25 min")


#############################################################################
#
# Test formatMinutesStr
#
#############################################################################


class test_getTimeInMinFromTotalTimeLine(unittest.TestCase):

  def test_None(self):
    self.assertEqual(
      getTimeInMinFromTotalTimeLine(
         "MPI_DEBUG", None),
      -1.0)

  def test_Empty(self):
    self.assertEqual(
      getTimeInMinFromTotalTimeLine(
         "MPI_DEBUG", ""),
      -1.0)

  def test_00(self):
    self.assertEqual(
      getTimeInMinFromTotalTimeLine(
         "MPI_DEBUG", "Total time for MPI_DEBUG = 1.16723643541 min"),
      1.16723643541)


#############################################################################
#
# Test extractPackageEnablesFromChangeStatus
#
#############################################################################


projectDepsXmlFileDefaultOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.gold.xml"
projectDependenciesDefault = getProjectDependenciesFromXmlFile(projectDepsXmlFileDefaultOverride)


class test_extractPackageEnablesFromChangeStatus(unittest.TestCase):


  def test_enable_all_and_other_packages(self):

    updateOutputStr = """
M	CMakeLists.txt
M	cmake/TrilinosPackages.cmake
M	cmake/python/checkin-test.py
M	doc/Thyra/coding_guildlines/ThyraCodingGuideLines.tex
P	packages/thyra/dummy.blah
A	packages/teuchos/example/ExplicitInstantiation/four_files/CMakeLists.txt
"""

    options = MockOptions()
    enablePackagesList = []

    extractPackageEnablesFromChangeStatus(updateOutputStr, options, GitRepo(""),
      enablePackagesList, False, projectDependenciesDefault)

    self.assertEqual( options.enableAllPackages, 'on' )
    self.assertEqual( enablePackagesList, [u'TrilinosFramework', u'Teuchos'] )


  def test_some_packages(self):

    updateOutputStr = """
? packages/triutils/doc/html
M	cmake/python/checkin-test.py
M	cmake/python/dump-cdash-deps-xml-file.py
A	packages/stratimikos/src/dummy.C
P       packages/stratimikos/dummy.blah
M	packages/thyra/src/Thyra_ConfigDefs.hpp
D	packages/tpetra/FSeconds.f
"""

    options = MockOptions()
    enablePackagesList = []

    extractPackageEnablesFromChangeStatus(updateOutputStr, options, GitRepo(""),
      enablePackagesList, False, projectDependenciesDefault)

    self.assertEqual( options.enableAllPackages, 'auto' )
    self.assertEqual( enablePackagesList, [u'TrilinosFramework', u'Stratimikos', u'ThyraCoreLibs', u'Tpetra'] )


  def test_extra_repo(self):

    updateOutputStr = """
M	CMakeLists.txt
M	ExtraTrilinosPackages.cmake
M	stalix/README
"""
    # NOTE: Above, we ignore top-level changes in extra repos which would cause global rebuilds
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    projectDependenciesLocal = getProjectDependenciesFromXmlFile(projectDepsXmlFileOverride)

    options = MockOptions()
    enablePackagesList = []

    extractPackageEnablesFromChangeStatus(updateOutputStr, options,
      GitRepo("preCopyrightTrilinos"),
      enablePackagesList, False, projectDependenciesLocal)

    self.assertEqual( options.enableAllPackages, 'auto' )
    self.assertEqual( enablePackagesList, [u'Stalix'] )



#############################################################################
#
#         Test getLastCommitMessageStrFromRawCommitLogStr
#
#############################################################################


class test_getLastCommitMessageStrFromRawCommitLogStr(unittest.TestCase):


  def test_clean_commit(self):
    cleanCommitMsg_expected = \
"""Some Commit Message

Some commit body

Some other message
"""
    rawLogOutput = "Standard git header stuff\n\n"+cleanCommitMsg_expected
    (cleanCommitMsg, numBlankLines) = getLastCommitMessageStrFromRawCommitLogStr(rawLogOutput)
    self.assertEqual(numBlankLines, -1)
    self.assertEqual(cleanCommitMsg, cleanCommitMsg_expected)


  def test_dirty_commit_1(self):
    cleanCommitMsg_expected = \
"""Some Commit Message

Some commit body

Some other message
"""
    rawLogOutput = \
       "Standard git header stuff\n\n" \
       +cleanCommitMsg_expected+ \
       "\nBuild/Test Cases Summary\n"
    #print "\nrawLogOutput:\n----------------\n", rawLogOutput, "----------------\n"
    (cleanCommitMsg, numBlankLines) = getLastCommitMessageStrFromRawCommitLogStr(rawLogOutput)
    self.assertEqual(numBlankLines, 1)
    self.assertEqual(cleanCommitMsg, cleanCommitMsg_expected)
    #print "\ncleanCommitMsg:\n----------------\n", cleanCommitMsg, "-----------------\n"


  def test_dirty_commit_2(self):
    cleanCommitMsg_expected = \
"""Some Commit Message

Some commit body

Some other message
"""
    rawLogOutput = \
       "Standard git header stuff\n\n" \
       +cleanCommitMsg_expected+ \
       "\nBuild/Test Cases Summary\n"
    #print "\nrawLogOutput:\n----------------\n", rawLogOutput, "----------------\n"
    (cleanCommitMsg, numBlankLines) = getLastCommitMessageStrFromRawCommitLogStr(rawLogOutput)
    self.assertEqual(numBlankLines, 1)
    self.assertEqual(cleanCommitMsg, cleanCommitMsg_expected)
    #print "\ncleanCommitMsg:\n----------------\n", cleanCommitMsg, "-----------------\n"


  def test_invalid_commit(self):
    cleanCommitMsg_expected = \
"""Some Commit Message

Some commit body

Some other message"""
    rawLogOutput = \
       "Standard git header stuff\n\n" \
       +cleanCommitMsg_expected+ \
       "\nBuild/Test Cases Summary\n"
    #print "\nrawLogOutput:\n----------------\n", rawLogOutput, "----------------\n"
    #(cleanCommitMsg, numBlankLines) = getLastCommitMessageStrFromRawCommitLogStr(rawLogOutput)
    self.assertRaises(Exception, getLastCommitMessageStrFromRawCommitLogStr, rawLogOutput)


  def test_two_summary_blocks(self):
    cleanCommitMsg_expected = \
"""Some Commit Message

Some commit body

Some other message
"""
    rawLogOutput = \
       "Standard git header stuff\n\n" \
       +cleanCommitMsg_expected+ \
       "\nBuild/Test Cases Summary\n"
    (cleanCommitMsg, numBlankLines) = getLastCommitMessageStrFromRawCommitLogStr(rawLogOutput)
    self.assertEqual(numBlankLines, 1)
    self.assertEqual(cleanCommitMsg, cleanCommitMsg_expected)
    # Strip it again to make sure we can pull it off again and recover
    rawLogOutput = \
       "Standard git header stuff\n\n" \
       +cleanCommitMsg+ \
       "\nBuild/Test Cases Summary\n"
    (cleanCommitMsg, numBlankLines) = getLastCommitMessageStrFromRawCommitLogStr(rawLogOutput)
    self.assertEqual(numBlankLines, 1)
    self.assertEqual(cleanCommitMsg, cleanCommitMsg_expected)


################################################################################
#
# Test Project name matching.
#
################################################################################

class test_matchProjectName(unittest.TestCase):
  def test_good_match(self):
    line = 'SET(PROJECT_NAME TestProject)'
    match = matchProjectName(line)
    self.assertEqual(match, 'TestProject')
    
  def test_match_with_extra_spaces(self):
    line = '  set ( PROJECT_NAME   TestProject ) '
    match = matchProjectName(line)
    self.assertEqual(match, 'TestProject')

  def test_no_match_wrong_variable(self):
    line = 'SET(SOME_VAR TestProject)'
    match = matchProjectName(line)
    self.assertFalse(match)

  def test_match_with_comment_at_end(self):
    line = 'Set(PROJECT_NAME TestProject) # This is a comment'
    match = matchProjectName(line)
    self.assertEqual(match, 'TestProject')


#############################################################################
#
# Test CMake helpers
#
#############################################################################

class test_cmakeScopedDefine(unittest.TestCase):
  def test_simple(self):
    result = cmakeScopedDefine('ProjectName', 'SOME_FLAG:BOOL', 'ON')
    self.assertEqual(result, '-DProjectName_SOME_FLAG:BOOL=ON')


#############################################################################
#
# Test TribitsGetExtraReposForCheckinTest.cmake 
#
#############################################################################


#run_extrarepo_test_verbose = True
run_extrarepo_test_verbose = False


def run_extrarepo_test(testObject, testName, extraReposFile, expectedReposList, \
  extraCmakeVars=None, expectedErrOutput=None \
  ):
  extraReposPythonOutFile = os.getcwd()+"/"+testName+".py"
  global g_withCmake
  cmnd = "\""+g_withCmake+"\""+ \
    " -DPROJECT_SOURCE_DIR="+mockProjectBaseDir+ \
    " -DTRIBITS_BASE_DIR="+tribitsBaseDir
  if extraCmakeVars:
    cmnd += " "+extraCmakeVars
  cmnd += \
    " -DEXTRA_REPOS_FILE="+scriptsDir+"/UnitTests/"+extraReposFile+ \
    " -DEXTRA_REPOS_PYTHON_OUT_FILE="+extraReposPythonOutFile+ \
    " -DUNITTEST_SKIP_FILTER_OR_ASSERT_EXTRA_REPOS=TRUE"+ \
    " -P "+scriptsDir+"/../package_arch/TribitsGetExtraReposForCheckinTest.cmake"
  consoleOutFile = testName+".out"
  rtn = echoRunSysCmnd(cmnd, throwExcept=False, timeCmnd=True, outFile=consoleOutFile,
    verbose=run_extrarepo_test_verbose)
  consoleOutputStr = readStrFromFile(consoleOutFile)
  if run_extrarepo_test_verbose:
    print "\nrtn =", rtn
    print "\n"+consoleOutFile+":\n", consoleOutputStr 
  if rtn == 0:
    readReposListTxt = readStrFromFile(extraReposPythonOutFile)
    if run_extrarepo_test_verbose:
      print "\nreadReposListTxt:\n", readReposListTxt
    readReposList = eval(readReposListTxt)
    if run_extrarepo_test_verbose:
      print "readReposList:\n", readReposList
    testObject.assertEqual(readReposList, expectedReposList)
  else:
    if run_extrarepo_test_verbose:
      print "\nexpectedErrOutput =", expectedErrOutput
    foundExpectedErrOutput = consoleOutputStr.find(expectedErrOutput)
    if foundExpectedErrOutput == -1:
      print "Error, failed to find:\n\n", expectedErrOutput
      print "\n\nin the output:\n\n", consoleOutputStr
    testObject.assertNotEqual(foundExpectedErrOutput, -1)


class test_TribitsGetExtraReposForCheckinTest(unittest.TestCase):

  def test_ExtraRepos1_implicit(self):
    run_extrarepo_test(
      self,
      "test_ExtraRepos1_implicit",
      "ExtraReposList_1.cmake",
      [
        {
          'NAME' : 'ExtraRepo1',
          'DIR' : 'ExtraRepo1',
          'REPOTYPE' : 'GIT',
          'REPOURL' : 'someurl.com:/git/data/SomeExtraRepo1',
          'PACKSTAT' : 'HASPACKAGES',
          'CATEGORY' : 'Continuous',
          },
        ],
      )

  def test_ExtraRepos1_explicit(self):
    run_extrarepo_test(
      self,
      "test_ExtraRepos1_explicit",
      "ExtraReposList_1.cmake",
      [
        {
          'NAME' : 'ExtraRepo1',
          'DIR' : 'ExtraRepo1',
          'REPOTYPE' : 'GIT',
          'REPOURL' : 'someurl.com:/git/data/SomeExtraRepo1',
          'PACKSTAT' : 'HASPACKAGES',
          'CATEGORY' : 'Continuous',
          },
        ],
      extraCmakeVars="-DENABLE_KNOWN_EXTERNAL_REPOS_TYPE=Continuous"
      )

  def test_ExtraReposAll(self):
    run_extrarepo_test(
      self,
      "test_ExtraReposAll",
      "ExtraReposList.cmake",
      [
        {
          'NAME' : 'ExtraRepo1',
          'DIR' : 'ExtraRepo1',
          'REPOTYPE' : 'GIT',
          'REPOURL' : 'someurl.com:/ExtraRepo1',
          'PACKSTAT' : 'HASPACKAGES',
          'CATEGORY' : 'Continuous',
          },
        {
          'NAME' : 'ExtraRepo2',
          'DIR' : 'packages/SomePackage/Blah',
          'REPOTYPE' : 'GIT',
          'REPOURL' : 'someurl2.com:/ExtraRepo2',
          'PACKSTAT' : 'NOPACKAGES',
          'CATEGORY' : 'Nightly',
          },
        {
          'NAME' : 'ExtraRepo3',
          'DIR' : 'ExtraRepo3',
          'REPOTYPE' : 'HG',
          'REPOURL' : 'someurl3.com:/ExtraRepo3',
          'PACKSTAT' : 'HASPACKAGES',
          'CATEGORY' : 'Continuous',
          },
        {
          'NAME' : 'ExtraRepo4',
          'DIR' : 'ExtraRepo4',
          'REPOTYPE' : 'SVN',
          'REPOURL' : 'someurl4.com:/ExtraRepo4',
          'PACKSTAT' : 'HASPACKAGES',
          'CATEGORY' : 'Nightly',
          },
        ],
      extraCmakeVars="-DENABLE_KNOWN_EXTERNAL_REPOS_TYPE=Nightly"
      )

  def test_ExtraReposContinuous(self):
    run_extrarepo_test(
      self,
      "test_ExtraReposContinuous",
      "ExtraReposList.cmake",
      [
        {
          'NAME' : 'ExtraRepo1',
          'DIR' : 'ExtraRepo1',
          'REPOTYPE' : 'GIT',
          'REPOURL' : 'someurl.com:/ExtraRepo1',
          'PACKSTAT' : 'HASPACKAGES',
          'CATEGORY' : 'Continuous',
          },
        {
          'NAME' : 'ExtraRepo3',
          'DIR' : 'ExtraRepo3',
          'REPOTYPE' : 'HG',
          'REPOURL' : 'someurl3.com:/ExtraRepo3',
          'PACKSTAT' : 'HASPACKAGES',
          'CATEGORY' : 'Continuous',
          },
        ],
      )

  def test_ExtraReposInvalidCategory(self):
    run_extrarepo_test(
      self,
      "test_ExtraReposInvalidCategory",
      "ExtraReposListInvalidCategory.cmake",
      [],
      )

  def test_ExtraReposInvalidType(self):
    run_extrarepo_test(
      self,
      "test_ExtraReposInvalidType",
      "ExtraReposListInvalidType.cmake",
      ["Will never get compared"],
      expectedErrOutput="Error, the repo type of 'InvalidType' for extra repo ExtraRepo1"
      )

  def test_ExtraReposEmptyList(self):
    run_extrarepo_test(
      self,
      "test_ExtraReposEmptyList",
      "ExtraReposListEmptyList.cmake",
      ["Will never get compared"],
      expectedErrOutput="Trilinos_EXTRAREPOS_DIR_REPOTYPE_REPOURL_PACKSTAT_CATEGORY is not defined!",
      )

  def test_ExtraReposEmptyFile(self):
    run_extrarepo_test(
      self,
      "test_ExtraReposEmptyFile",
      "ExtraReposListEmptyFile.cmake",
      ["Will never get compared"],
      expectedErrOutput="Trilinos_EXTRAREPOS_DIR_REPOTYPE_REPOURL_PACKSTAT_CATEGORY is not defined!",
      )


#############################################################################
#
# Test TribitsGitRepos
#
#############################################################################


def assertCompareGitRepoLists(testObject, gitRepoList, expectedGitRepoList):
  testObject.assertEqual(len(gitRepoList), len(expectedGitRepoList))
  for i in range(len(gitRepoList)):
    gitRepo = gitRepoList[i]
    expectedGitRepo = expectedGitRepoList[i]
    testObject.assertEqual(gitRepo.repoName, expectedGitRepo.repoName)
    testObject.assertEqual(gitRepo.repoDir, expectedGitRepo.repoDir)
    testObject.assertEqual(gitRepo.repoType, expectedGitRepo.repoType)
    testObject.assertEqual(gitRepo.repoHasPackages, expectedGitRepo.repoHasPackages)
    testObject.assertEqual(gitRepo.hasChanges, expectedGitRepo.hasChanges)


def assertCompareTribitGitRepos(testObject, tribitsGitRepo, expectedTribitsGitRepo):
  assertCompareGitRepoLists(testObject, tribitsGitRepo.gitRepoList(),
     expectedTribitsGitRepo.gitRepoList())
  testObject.assertEqual(tribitsGitRepo.tribitsExtraRepoNamesList(),
    expectedTribitsGitRepo.tribitsExtraRepoNamesList())


#test_TribitsGitRepos_verbose = True
test_TribitsGitRepos_verbose = False


def test_TribitsGitRepos_run_case(testObject, testName, inOptions, \
  expectPass, \
  expectedTribitsExtraRepoNamesList, expectedGitRepos, \
  consoleRegexMatches=None, consoleRegexNotMatches=None, \
  exceptionRegexMatches=None \
  ):
  inOptions.withCmake = g_withCmake
  currDir = os.getcwd()
  if os.path.exists(testName):
    runSysCmnd("rm -rf "+testName)
  os.mkdir(testName)
  os.chdir(testName)
  try:
    consoleOutputFile = "Console.out" 
    tribitsGitRepos = TribitsGitRepos()
    cmndPassed = False
    try:
      tribitsGitRepos.initFromCommandlineArguments(inOptions, \
        consoleOutputFile=consoleOutputFile, \
        verbose=test_TribitsGitRepos_verbose)
      cmndPassed = True
      # NOTE: the file consoleOutputFile still gets written, even if throw
    except Exception, e:
      #print "e =", e
      if exceptionRegexMatches:
        eMsg = e.args[0]
        for exceptRegex in exceptionRegexMatches.split('\n'):
          matchResult = re.search(exceptRegex, eMsg)
          if not matchResult:
            print "Error, the regex expression '"+exceptRegex+"' was not" \
              +" found in the exception string '"+eMsg+"'!"
          testObject.assertNotEqual(matchResult, None)
    testObject.assertEqual(cmndPassed, expectPass)
    if cmndPassed:
      #print "\ntribitsGitRepos =", tribitsGitRepos
      testObject.assertEqual(tribitsGitRepos.numTribitsExtraRepos(), len(expectedTribitsExtraRepoNamesList))
      expectedTribitsGitRepo = TribitsGitRepos().reset()
      expectedTribitsGitRepo._TribitsGitRepos__gitRepoList.extend(expectedGitRepos)
      testObject.assertEqual(tribitsGitRepos.tribitsExtraRepoNamesList(), expectedTribitsExtraRepoNamesList)
      expectedTribitsGitRepo._TribitsGitRepos__tribitsExtraRepoNamesList.extend(expectedTribitsExtraRepoNamesList)
      assertCompareTribitGitRepos(testObject, tribitsGitRepos, expectedTribitsGitRepo)
    if consoleRegexMatches:
      assertGrepFileForRegexStrList(testObject, testName, consoleOutputFile,
        consoleRegexMatches, test_TribitsGitRepos_verbose)
    if consoleRegexNotMatches:
      assertNotGrepFileForRegexStrList(testObject, testName, consoleOutputFile,
        consoleRegexNotMatches, test_TribitsGitRepos_verbose)
  finally:
    os.chdir(currDir)


class test_TribitsGitRepos(unittest.TestCase):

  def test_noExtraRepos(self):
    tribitsGitRepos = TribitsGitRepos()
    expectedTribitsGitRepo = TribitsGitRepos().reset()
    expectedTribitsGitRepo._TribitsGitRepos__gitRepoList.append(GitRepo("", "", "GIT", True))
    self.assertEqual(tribitsGitRepos.tribitsExtraRepoNamesList(), [])
    self.assertEqual(tribitsGitRepos.numTribitsExtraRepos(), 0)
    assertCompareTribitGitRepos(self, tribitsGitRepos, expectedTribitsGitRepo)

  def test_noExtraReposFile_extraRepos(self):
    testName = "test_noExtraReposFile_extraRepos"
    inOptions = MockOptions()
    inOptions.extraRepos = "preCopyrightTrilinos"
    inOptions.extraReposFile = ""
    expectedTribitsExtraRepoNamesList = ["preCopyrightTrilinos"]
    expectedGitRepos = [
      GitRepo("", "", "GIT", True),
      GitRepo('preCopyrightTrilinos', 'preCopyrightTrilinos', "GIT", True),
      ]
    consoleRegexMatches = None
    consoleRegexNotMatches = None
    test_TribitsGitRepos_run_case(self, testName, inOptions, True, \
      expectedTribitsExtraRepoNamesList, expectedGitRepos, \
      consoleRegexMatches, consoleRegexNotMatches)

  def test_ExtraRepos3_Continuous(self):
    testName = "test_ExtraRepos3_Continuous"
    inOptions = MockOptions()
    inOptions.extraReposFile = \
      tribitsBaseDir+"/python/UnitTests/ExtraReposListExisting_3.cmake"
    inOptions.extraReposType = "Continuous"
    expectedTribitsExtraRepoNamesList = ["preCopyrightTrilinos"]
    expectedGitRepos = [
      GitRepo("", "", "GIT", True),
      GitRepo('preCopyrightTrilinos', 'preCopyrightTrilinos', "GIT", True),
      GitRepo('ExtraTeuchosRepo', 'packages/teuchos/extrastuff', 'GIT', False)
      ]
    consoleRegexMatches = \
      "Adding extra Continuous repository preCopyrightTrilinos\n"+\
      "Adding extra Continuous repository ExtraTeuchosRepo\n"
    consoleRegexNotMatches = \
      "Adding extra Nightly repository extraTrilinosRepo\n"
    test_TribitsGitRepos_run_case(self, testName, inOptions, True, \
      expectedTribitsExtraRepoNamesList, expectedGitRepos, \
      consoleRegexMatches, consoleRegexNotMatches)

  def test_ExtraRepos3_Nightly(self):
    testName = "test_ExtraRepos3_Nightly"
    inOptions = MockOptions()
    inOptions.extraReposFile = \
      tribitsBaseDir+"/python/UnitTests/ExtraReposListExisting_3.cmake"
    inOptions.extraReposType = "Nightly"
    expectedPass = True
    expectedTribitsExtraRepoNamesList = ["preCopyrightTrilinos", "extraTrilinosRepo"]
    expectedGitRepos = [
      GitRepo("", "", "GIT", True),
      GitRepo('preCopyrightTrilinos', 'preCopyrightTrilinos', "GIT", True),
      GitRepo('extraTrilinosRepo', 'extraTrilinosRepo', "GIT", True),
      GitRepo('ExtraTeuchosRepo', 'packages/teuchos/extrastuff', 'GIT', False)
      ]
    consoleRegexMatches = \
      "Adding extra Continuous repository preCopyrightTrilinos\n"+\
      "Adding extra Continuous repository ExtraTeuchosRepo\n"+\
      "Adding extra Nightly repository extraTrilinosRepo\n"
    consoleRegexNotMatches = None
    test_TribitsGitRepos_run_case(self, testName, inOptions, expectedPass, \
      expectedTribitsExtraRepoNamesList, expectedGitRepos, \
      consoleRegexMatches, consoleRegexNotMatches)

  def test_ExtraReposExisting1Missing1_assert(self):
    testName = "test_ExtraReposExisting1Missing1_assert"
    inOptions = MockOptions()
    inOptions.extraReposFile = \
      tribitsBaseDir+"/python/UnitTests/ExtraReposListExisting1Missing1.cmake"
    inOptions.extraReposType = "Nightly"
    expectedPass = False
    expectedTribitsExtraRepoNamesList = ["Will never be compared"]
    expectedGitRepos = ["Will never be compared."]
    consoleRegexMatches = \
      "ERROR! Skipping missing extra repo .MissingRepo. since\n"
    consoleRegexNotMatches = \
      "Adding extra Continuous repository MissingRepo\n"
    test_TribitsGitRepos_run_case(self, testName, inOptions, expectedPass, \
      expectedTribitsExtraRepoNamesList, expectedGitRepos, \
    consoleRegexMatches, consoleRegexNotMatches)

  def test_ExtraReposExisting1Missing1_ignore(self):
    testName = "test_ExtraReposExisting1Missing1_ignore"
    inOptions = MockOptions()
    inOptions.extraReposFile = \
      tribitsBaseDir+"/python/UnitTests/ExtraReposListExisting1Missing1.cmake"
    inOptions.extraReposType = "Nightly"
    inOptions.ignoreMissingExtraRepos = True
    expectedPass = True
    expectedTribitsExtraRepoNamesList = ["preCopyrightTrilinos"]
    expectedGitRepos = [
      GitRepo("", "", "GIT", True),
      GitRepo('preCopyrightTrilinos', 'preCopyrightTrilinos', "GIT", True),
      ]
    consoleRegexMatches = \
      "WARNING!  Ignoring missing extra repo .MissingRepo. as requested since\n"
    consoleRegexNotMatches = \
      "Adding extra Continuous repository MissingRepo"
    test_TribitsGitRepos_run_case(self, testName, inOptions, expectedPass, \
      expectedTribitsExtraRepoNamesList, expectedGitRepos, \
    consoleRegexMatches, consoleRegexNotMatches)

  def test_ExtraRepos3_listExtraReposNotListed1(self):
    testName = "test_ExtraRepos3_listExtraReposNotListed1"
    inOptions = MockOptions()
    inOptions.extraReposFile = \
      tribitsBaseDir+"/python/UnitTests/ExtraReposListExisting_3.cmake"
    inOptions.extraRepos = "extraRepoNotInList"
    inOptions.extraReposType = "Nightly"
    expectedPass = False
    expectedTribitsExtraRepoNamesList = ["Will never be compared"]
    expectedGitRepos = ["Will never be compared"]
    consoleRegexMatches = \
      "ERROR! The list of extra repos passed in .extraRepoNotInList. is not\n"
    consoleRegexNotMatches = None
    test_TribitsGitRepos_run_case(self, testName, inOptions, expectedPass, \
      expectedTribitsExtraRepoNamesList, expectedGitRepos, \
      consoleRegexMatches, consoleRegexNotMatches)

  def test_ExtraRepos3_listExtraReposNotListed2(self):
    testName = "test_ExtraRepos3_listExtraReposNotListed2"
    inOptions = MockOptions()
    inOptions.extraReposFile = \
      tribitsBaseDir+"/python/UnitTests/ExtraReposListExisting_3.cmake"
    inOptions.extraRepos = "preCopyrightTrilinos,extraRepoNotInList"
    inOptions.extraReposType = "Nightly"
    expectedPass = False
    expectedTribitsExtraRepoNamesList = ["Will never be compared"]
    expectedGitRepos = ["Will never be compared"]
    consoleRegexMatches = \
      "ERROR! The list of extra repos passed in\n"+\
     ".preCopyrightTrilinos.extraRepoNotInList. is not a subset and in the same\n"
    consoleRegexNotMatches = None
    test_TribitsGitRepos_run_case(self, testName, inOptions, expectedPass, \
      expectedTribitsExtraRepoNamesList, expectedGitRepos, \
      consoleRegexMatches, consoleRegexNotMatches)

  def test_ExtraRepos3_extraReposFullList_right_order(self):
    testName = "test_ExtraRepos3_extraReposFullList_right_order"
    inOptions = MockOptions()
    inOptions.extraReposFile = \
      tribitsBaseDir+"/python/UnitTests/ExtraReposListExisting_3.cmake"
    inOptions.extraReposType = "Nightly"
    inOptions.extraRepos = "preCopyrightTrilinos,extraTrilinosRepo,ExtraTeuchosRepo"
    expectedPass = True
    expectedTribitsExtraRepoNamesList = ["preCopyrightTrilinos", "extraTrilinosRepo"]
    expectedGitRepos = [
      GitRepo("", "", "GIT", True),
      GitRepo('preCopyrightTrilinos', 'preCopyrightTrilinos', "GIT", True),
      GitRepo('extraTrilinosRepo', 'extraTrilinosRepo', "GIT", True),
      GitRepo('ExtraTeuchosRepo', 'packages/teuchos/extrastuff', 'GIT', False)
      ]
    consoleRegexMatches = None
    consoleRegexNotMatches = None
    test_TribitsGitRepos_run_case(self, testName, inOptions, expectedPass, \
      expectedTribitsExtraRepoNamesList, expectedGitRepos, \
      consoleRegexMatches, consoleRegexNotMatches)

  def test_ExtraRepos3_extraReposFullList_wrong_order(self):
    testName = "test_ExtraRepos3_extraReposFullList_wrong_order"
    inOptions = MockOptions()
    inOptions.extraReposFile = \
      tribitsBaseDir+"/python/UnitTests/ExtraReposListExisting_3.cmake"
    inOptions.extraReposType = "Nightly"
    inOptions.extraRepos = "extraTrilinosRepo,preCopyrightTrilinos"
    expectedPass = False
    expectedTribitsExtraRepoNamesList = ["Will never be compared"]
    expectedGitRepos = ["Will never be compared"]
    consoleRegexMatches = \
      "ERROR! The list of extra repos passed in\n"+\
     ".extraTrilinosRepo;preCopyrightTrilinos. is not a subset and in the same\n"
    consoleRegexNotMatches = None
    test_TribitsGitRepos_run_case(self, testName, inOptions, expectedPass, \
      expectedTribitsExtraRepoNamesList, expectedGitRepos, \
      consoleRegexMatches, consoleRegexNotMatches)

  def test_ExtraRepos3_listExtraRepos1_first(self):
    testName = "test_ExtraRepos3_listExtraRepos1_first"
    inOptions = MockOptions()
    inOptions.extraReposFile = \
      tribitsBaseDir+"/python/UnitTests/ExtraReposListExisting_3.cmake"
    inOptions.extraRepos = "preCopyrightTrilinos"
    inOptions.extraReposType = "Nightly"
    expectedPass = True
    expectedTribitsExtraRepoNamesList = ["preCopyrightTrilinos"]
    expectedGitRepos = [
      GitRepo("", "", "GIT", True),
      GitRepo('preCopyrightTrilinos', 'preCopyrightTrilinos', "GIT", True),
      ]
    consoleRegexMatches = None
    consoleRegexNotMatches = None
    test_TribitsGitRepos_run_case(self, testName, inOptions, expectedPass, \
      expectedTribitsExtraRepoNamesList, expectedGitRepos, \
      consoleRegexMatches, consoleRegexNotMatches)

  def test_ExtraRepos3_listExtraRepos1_middle(self):
    testName = "test_ExtraRepos3_listExtraRepos1_middle"
    inOptions = MockOptions()
    inOptions.extraReposFile = \
      tribitsBaseDir+"/python/UnitTests/ExtraReposListExisting_3.cmake"
    inOptions.extraRepos = "extraTrilinosRepo"
    inOptions.extraReposType = "Nightly"
    expectedPass = True
    expectedTribitsExtraRepoNamesList = ["extraTrilinosRepo"]
    expectedGitRepos = [
      GitRepo("", "", "GIT", True),
      GitRepo('extraTrilinosRepo', 'extraTrilinosRepo', "GIT", True),
      ]
    consoleRegexMatches = None
    consoleRegexNotMatches = None
    test_TribitsGitRepos_run_case(self, testName, inOptions, expectedPass, \
      expectedTribitsExtraRepoNamesList, expectedGitRepos, \
      consoleRegexMatches, consoleRegexNotMatches)

  def test_ExtraRepos3_listExtraRepos1_last(self):
    testName = "test_ExtraRepos3_listExtraRepos1_last"
    inOptions = MockOptions()
    inOptions.extraReposFile = \
      tribitsBaseDir+"/python/UnitTests/ExtraReposListExisting_3.cmake"
    inOptions.extraRepos = "ExtraTeuchosRepo"
    inOptions.extraReposType = "Nightly"
    expectedPass = True
    expectedTribitsExtraRepoNamesList = []
    expectedGitRepos = [
      GitRepo("", "", "GIT", True),
      GitRepo('ExtraTeuchosRepo', 'packages/teuchos/extrastuff', 'GIT', False)
      ]
    consoleRegexMatches = None
    consoleRegexNotMatches = None
    test_TribitsGitRepos_run_case(self, testName, inOptions, expectedPass, \
      expectedTribitsExtraRepoNamesList, expectedGitRepos, \
      consoleRegexMatches, consoleRegexNotMatches)

  def test_ExtraRepos3NoContinuous_noExtraRepos(self):
    testName = "test_ExtraRepos3NoContinuous_noExtraRepos"
    inOptions = MockOptions()
    inOptions.extraReposFile = \
      tribitsBaseDir+"/python/UnitTests/ExtraReposList3NoContinuous.cmake"
    inOptions.extraRepos = ""
    inOptions.extraReposType = "Continuous"
    expectedPass = True
    expectedTribitsExtraRepoNamesList = []
    expectedGitRepos = [
      GitRepo("", "", "GIT", True),
      ]
    consoleRegexMatches = None
    consoleRegexNotMatches = None
    test_TribitsGitRepos_run_case(self, testName, inOptions, expectedPass, \
      expectedTribitsExtraRepoNamesList, expectedGitRepos, \
      consoleRegexMatches, consoleRegexNotMatches)

  def test_ExtraReposHasPackagesAndDeepDir(self):
    testName = "test_ExtraReposHasPackagesAndDeepDir"
    inOptions = MockOptions()
    inOptions.extraReposFile = \
      tribitsBaseDir+"/python/UnitTests/ExtraReposListHasPackagesAndDeepDir.cmake"
    inOptions.extraRepos = ""
    inOptions.extraReposType = "Continuous"
    expectedPass = False
    expectedTribitsExtraRepoNamesList = [ "Will never be compared"]
    expectedGitRepos = ["Will never be compared"]
    consoleRegexMatches = None
    consoleRegexNotMatches = None
    exceptionRegexMatches = \
      "ERROR!  For extra repo 'ExtraTeuchosRepo', if repoHasPackages==True then repoDir must be same as repo name, not 'packages/teuchos/extrastuff'!\n"
    test_TribitsGitRepos_run_case(self, testName, inOptions, expectedPass, \
      expectedTribitsExtraRepoNamesList, expectedGitRepos, \
      consoleRegexMatches, consoleRegexNotMatches, exceptionRegexMatches)

  def test_ExtraReposNotGit(self):
    testName = "test_ExtraReposNotGit"
    inOptions = MockOptions()
    inOptions.extraReposFile = \
      tribitsBaseDir+"/python/UnitTests/ExtraReposListNotGit.cmake"
    inOptions.extraRepos = ""
    inOptions.extraReposType = "Continuous"
    expectedPass = False
    expectedTribitsExtraRepoNamesList = [ "Will never be compared"]
    expectedGitRepos = ["Will never be compared"]
    consoleRegexMatches = None
    consoleRegexNotMatches = None
    exceptionRegexMatches = \
      "ERROR!  For extra repo 'ExtraTeuchosRepo', the repo type 'SVN' is not supported by the checkin-test.py script, only 'GIT'!\n"
    test_TribitsGitRepos_run_case(self, testName, inOptions, expectedPass, \
      expectedTribitsExtraRepoNamesList, expectedGitRepos, \
      consoleRegexMatches, consoleRegexNotMatches, exceptionRegexMatches)

  def test_ExraRepoListProjectDefault(self):
    testName = "test_ExraRepoListProjectDefault"
    inOptions = MockOptions()
    inOptions.extraReposFile = "project"
    inOptions.extraRepos = "preCopyrightTrilinos,extraTrilinosRepo"
    inOptions.extraReposType = "Nightly"
    expectedPass = True
    expectedTribitsExtraRepoNamesList = ["preCopyrightTrilinos", "extraTrilinosRepo"]
    expectedGitRepos = [
      GitRepo("", "", "GIT", True),
      GitRepo('preCopyrightTrilinos', 'preCopyrightTrilinos', "GIT", True),
      GitRepo('extraTrilinosRepo', 'extraTrilinosRepo', "GIT", True),
      ]
    consoleRegexMatches = \
      "Adding extra Continuous repository preCopyrightTrilinos\n"+\
      "Adding extra Nightly repository extraTrilinosRepo\n"
    consoleRegexNotMatches = None
    exceptionRegexMatches = None
    test_TribitsGitRepos_run_case(self, testName, inOptions, expectedPass, \
      expectedTribitsExtraRepoNamesList, expectedGitRepos, \
      consoleRegexMatches, consoleRegexNotMatches, exceptionRegexMatches)



#############################################################################
#
# Test checkin-test.py script
#
#############################################################################


# Test Data

g_cmndinterceptsCurrentBranch = \
  "IT: eg branch; 0; '* currentbranch'\n"

g_cmndinterceptsStatusPasses = \
  "IT: eg status; 0; '(on master branch)'\n"

g_cmndinterceptsStatusChangedButNotUpdatedPasses = \
  "IT: eg status; 0; 'Changed but not updated'\n"

g_cmndinterceptsPullOnlyPasses = \
  "IT: eg pull; 0; 'pulled changes passes'\n"

g_cmndinterceptsPullOnlyFails = \
  "IT: eg pull; 1; 'pull failed'\n"

g_cmndinterceptsPullOnlyNoUpdatesPasses = \
  "IT: eg pull; 0; 'Already up-to-date.'\n"

g_cmndinterceptsStatusPullPasses = \
  g_cmndinterceptsStatusPasses+ \
  g_cmndinterceptsPullOnlyPasses

g_cmndinterceptsDiffOnlyPasses = \
  "IT: eg diff --name-status origin/currentbranch; 0; 'M\tpackages/teuchos/CMakeLists.txt'\n"

g_cmndinterceptsDiffOnlyNoChangesPasses = \
  "IT: eg diff --name-status origin/currentbranch; 0; ''\n"

g_cmndinterceptsDiffOnlyPassesPreCopyrightTrilinos = \
  "IT: eg diff --name-status origin/currentbranch; 0; 'M\tteko/CMakeLists.txt'\n"

g_cmndinterceptsDiffOnlyPassesExtraTrilinosRepo = \
  "IT: eg diff --name-status origin/currentbranch; 0; 'M\textrapack/src/ExtraPack_ConfigDefs.hpp'\n"

g_cmndinterceptsPullPasses = \
  g_cmndinterceptsStatusPullPasses \
  +g_cmndinterceptsDiffOnlyPasses

g_cmndinterceptsConfigPasses = \
  "IT: \./do-configure; 0; 'do-configure passed'\n"

g_cmndinterceptsConfigBuildPasses = \
  g_cmndinterceptsConfigPasses+ \
  "IT: make -j3; 0; 'make passed'\n"

g_cmndinterceptsConfigBuildTestPasses = \
  g_cmndinterceptsConfigBuildPasses+ \
  "IT: ctest -j5; 0; '100% tests passed, 0 tests failed out of 100'\n"

g_cmnginterceptsEgLogCmnds = \
  "IT: eg cat-file -p HEAD; 0; 'This is the last commit message'\n" \
  "IT: eg log --oneline currentbranch \^origin/currentbranch; 0; '12345 Only one commit'\n" \
  "IT: eg log --pretty=format:'%h' currentbranch\^ \^origin/currentbranch; 0; '12345'\n"

g_cmndinterceptsFinalPullRebasePasses = \
  "IT: eg pull && eg rebase --against origin/currentbranch; 0; 'final eg pull and rebase passed'\n"

g_cmndinterceptsFinalPullRebaseFails = \
  "IT: eg pull && eg rebase --against origin/currentbranch; 1; 'final eg pull and rebase failed'\n"

g_cmndinterceptsAmendCommitPasses = \
  "IT: eg commit --amend -F .*; 0; 'Amending the last commit passed'\n"

g_cmndinterceptsAmendCommitFails = \
  "IT: eg commit --amend -F .*; 1; 'Amending the last commit failed'\n"

g_cmndinterceptsLogCommitsPasses = \
  "IT: eg log --oneline currentbranch \^origin/currentbranch; 0; '54321 Only one commit'\n"

g_cmndinterceptsPushOnlyPasses = \
  "IT: eg push; 0; 'push passes'\n"

g_cmndinterceptsPushOnlyFails = \
  "IT: eg push; 1; 'push failed'\n"

g_cmndinterceptsFinalPushPasses = \
  g_cmndinterceptsFinalPullRebasePasses+\
  g_cmnginterceptsEgLogCmnds+ \
  g_cmndinterceptsAmendCommitPasses+ \
  g_cmndinterceptsLogCommitsPasses+ \
  "IT: cat modifiedFiles.out; 0; 'M\tpackages/teuchos/CMakeLists.txt'\n"\
  "IT: eg push; 0; 'push passes'\n"

g_cmndinterceptsCatModifiedFilesPasses = \
  "IT: cat modifiedFiles.out; 0; 'M\tpackages/teuchos/CMakeLists.txt'\n"

g_cmndinterceptsCatModifiedFilesNoChanges = \
  "IT: cat modifiedFiles.out; 0; ''\n"

g_cmndinterceptsCatModifiedFilesPreCoprightTrilinosPasses = \
  "IT: cat modifiedFiles.preCopyrightTrilinos.out; 0; 'M\tteko/CMakeLists.txt'\n"

g_cmndinterceptsCatModifiedFilesPreCoprightTrilinosNoChanges = \
  "IT: cat modifiedFiles.preCopyrightTrilinos.out; 0; ''\n"

g_cmndinterceptsFinalPushNoAppendTestResultsPasses = \
  "IT: eg pull && eg rebase --against origin/currentbranch; 0; 'final eg pull and rebase passed'\n" \
  +g_cmndinterceptsLogCommitsPasses\
  +g_cmndinterceptsCatModifiedFilesPasses\
  +g_cmndinterceptsPushOnlyPasses

g_cmndinterceptsFinalPushNoRebasePasses = \
  "IT: eg pull; 0; 'final eg pull only passed'\n" \
  +g_cmnginterceptsEgLogCmnds+ \
  "IT: eg commit --amend -F .*; 0; 'Amending the last commit passed'\n" \
  +g_cmndinterceptsLogCommitsPasses\
  +g_cmndinterceptsCatModifiedFilesPasses\
  +g_cmndinterceptsPushOnlyPasses

g_cmndinterceptsSendBuildTestCaseEmail = \
  "IT: mailx -s .*; 0; 'Do not really sending build/test case email'\n"

g_cmndinterceptsSendFinalEmail = \
  "IT: mailx -s .*; 0; 'Do not really send email '\n"

g_cmndinterceptsExtraRepo1ThroughStatusPasses = \
  "IT: .*cmake .+ -P .+/TribitsDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
  +g_cmndinterceptsCurrentBranch \
  +g_cmndinterceptsStatusPasses \
  +g_cmndinterceptsStatusPasses

g_cmndinterceptsExtraRepo1DoAllThroughTest = \
  g_cmndinterceptsExtraRepo1ThroughStatusPasses \
  +g_cmndinterceptsPullOnlyPasses \
  +g_cmndinterceptsPullOnlyPasses \
  +g_cmndinterceptsDiffOnlyPasses \
  +g_cmndinterceptsDiffOnlyPassesPreCopyrightTrilinos \
  +g_cmndinterceptsConfigBuildTestPasses \
  +g_cmndinterceptsSendBuildTestCaseEmail

g_cmndinterceptsExtraRepo1DoAllUpToPush = \
  g_cmndinterceptsExtraRepo1DoAllThroughTest \
  +g_cmndinterceptsFinalPullRebasePasses \
  +g_cmndinterceptsFinalPullRebasePasses \
  +g_cmnginterceptsEgLogCmnds \
  +g_cmndinterceptsAmendCommitPasses \
  +g_cmnginterceptsEgLogCmnds \
  +g_cmndinterceptsAmendCommitPasses \
  +g_cmndinterceptsLogCommitsPasses \
  +g_cmndinterceptsLogCommitsPasses

g_expectedRegexUpdatePasses = \
  "Update passed!\n" \

g_expectedRegexUpdateWithBuildCasePasses = \
  "Update passed!\n" \
  "The update passed!\n" \
  "Update: Passed\n"

g_expectedRegexConfigPasses = \
  "Modified file: .packages/teuchos/CMakeLists\.txt\n" \
  "  => Enabling .Teuchos.!\n" \
  "Configure passed!\n" \
  "The configure passed!\n" \
  "Configure: Passed\n" \

g_expectedRegexBuildPasses = \
  "Build passed!\n" \
  "The build passed!\n" \
  "Build: Passed\n"

g_expectedRegexBuildFailed = \
  "Build failed returning 1!\n" \
  "The build FAILED!\n" \
  "Build: FAILED\n"

g_expectedRegexTestPasses = \
  "No tests failed!\n" \
  "testResultsLine = .100% tests passed, 0 tests failed out of 100.\n" \
  "passed: passed=100,notpassed=0\n" \
  "Test: Passed\n"

g_expectedRegexTestNotRun = \
  "The tests where never even run!\n" \
  "Test: FAILED\n"

g_expectedCommonOptionsSummary = \
  "Enabled Packages: Teuchos\n" \
  "Make Options: -j3\n" \
  "CTest Options: -j5\n"

g_verbose=True
g_verbose=False


#
# Test helper functions
#


g_checkin_test_tests_dir = "checkin_test_tests"


def create_checkin_test_case_dir(testName, verbose=False):
  baseDir = os.getcwd()
  testDirName = os.path.join(g_checkin_test_tests_dir, testName)
  createDir(g_checkin_test_tests_dir, verbose)
  createDir(testDirName, verbose)
  return testDirName


# Main unit test driver
def checkin_test_run_case(testObject, testName, optionsStr, cmndInterceptsStr, \
  expectPass, passRegexStrList, filePassRegexStrList=None, mustHaveCheckinTestOut=True, \
  failRegexStrList=None, fileFailRegexStrList=None, envVars=[], inPathEg=True, egVersion=True \
  ):

  verbose = g_verbose

  passRegexList = passRegexStrList.split('\n')

  if verbose: print "\npassRegexList =", passRegexList

  # A) Create the test directory

  baseDir = os.getcwd()
  echoChDir(create_checkin_test_case_dir(testName, verbose), verbose)

  try:

    # B) Create the command to run the checkin-test.py script

    
    cmndArgs = [
      scriptsDir + "/../checkin-test.py",
      "--with-cmake=\""+g_withCmake+"\"",
      "--project-name=Trilinos",
      "--no-eg-git-version-check",
      "--src-dir="+scriptsDir+"/../package_arch/UnitTests/MockTrilinos",
      "--send-email-to=bogous@somwhere.com",
      "--project-configuration=%s" % os.path.join(scriptsDir,
        'UnitTests', 'CheckinTest_UnitTests_Config.py'),
      optionsStr,
      ]
    cmnd = ' '.join(cmndArgs)
    # NOTE: Above, we want to turn off the eg/git version tests since we want
    # these unit tests to run on machines that do not have the official
    # versions (e.g. the SCICO LAN) but where the versions might be okay.
    # Also, we have to point to the static mock Trilinos source directory
    # also so that preCopyrighTrilinos will show up as an extra repo.
    
    # C) Set up the command intercept file

    baseCmndInterceptsStr = \
      "FT: .*checkin-test-impl\.py.*\n" \
      "FT: .*cmake .*TribitsGetExtraReposForCheckinTest.cmake.*\n" \
      "FT: date\n" \
      "FT: rm [a-zA-Z0-9_/\.]+\n" \
      "FT: touch .*\n" \
      "FT: chmod .*\n" \
      "FT: hostname\n" \
      "FT: grep .*"+getTestOutputFileName()+"\n" \
      "FT: grep .*"+getEmailBodyFileName()+"\n" \
      "FT: grep .*REQUESTED ACTIONS\: PASSED.*\n"

    if inPathEg:
      baseCmndInterceptsStr += \
      "IT: git config --get user.email; 0; bogous@somwhere.com\n" \
      +"IT: which eg; 0; /some/path/eg\n"

    if egVersion:
      baseCmndInterceptsStr += \
      "IT: eg --version; 0; "+g_officialEgVersion+"\n"

    fullCmndInterceptsStr = baseCmndInterceptsStr + cmndInterceptsStr

    fullCmndInterceptsFileName = os.path.join(os.getcwd(), "cmndIntercepts.txt")
    writeStrToFile(fullCmndInterceptsFileName, fullCmndInterceptsStr)

    os.environ['GENERAL_SCRIPT_SUPPORT_CMND_INTERCEPTS_FILE'] = fullCmndInterceptsFileName

    os.environ['CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE'] = projectDepsXmlFileDefaultOverride
    
    # D) Run the checkin-test.py script with mock commands

    for envVar in envVars:
      (varName,varValue) = envVar.split("=")
      #print "varName="+varName
      #print "varValue="+varValue
      os.environ[varName] = varValue
      
    checkin_test_test_out = "checkin-test.test.out"

    rtnCode = echoRunSysCmnd(cmnd, timeCmnd=True, throwExcept=False,
      outFile=checkin_test_test_out, verbose=verbose)
    
    # E) Grep the main output file looking for specific strings

    if mustHaveCheckinTestOut:
      outputFileToGrep = "checkin-test.out"
    else:
      outputFileToGrep = checkin_test_test_out

    assertGrepFileForRegexStrList(testObject, testName, outputFileToGrep,
      passRegexStrList, verbose)

    if failRegexStrList:
      assertNotGrepFileForRegexStrList(testObject, testName, outputFileToGrep,
        failRegexStrList, verbose)

    # F) Grep a set of output files looking for given strings

    if filePassRegexStrList:
      for fileRegexGroup in filePassRegexStrList:
        (fileName, regexStrList) = fileRegexGroup
        assertGrepFileForRegexStrList(testObject, testName, fileName, regexStrList, verbose)

    if fileFailRegexStrList:
      for fileRegexGroup in fileFailRegexStrList:
        (fileName, regexStrList) = fileRegexGroup
        assertNotGrepFileForRegexStrList(testObject, testName, fileName, regexStrList, verbose)

    # G) Examine the final return code

    if expectPass:
      testObject.assertEqual(rtnCode, 0)
    else:
      testObject.assertNotEqual(rtnCode, 0)
    
  finally:
    # \H) Get back to the current directory and reset
    echoChDir(baseDir, verbose=verbose)
    os.environ['GENERAL_SCRIPT_SUPPORT_CMND_INTERCEPTS_FILE']=""


# Helper test case that is used as the inital case for other tests
def g_test_do_all_default_builds_mpi_debug_pass(testObject, testName):
  checkin_test_run_case(
    \
    testObject,
    \
    testName,
    \
    "--make-options=-j3 --ctest-options=-j5 --default-builds=MPI_DEBUG --do-all",
    \
    g_cmndinterceptsCurrentBranch \
    +g_cmndinterceptsPullPasses \
    +g_cmndinterceptsConfigBuildTestPasses \
    +g_cmndinterceptsSendBuildTestCaseEmail \
    +g_cmndinterceptsSendFinalEmail \
    ,
    \
    True,
    \
    g_expectedRegexUpdateWithBuildCasePasses \
    +g_expectedRegexConfigPasses \
    +g_expectedRegexBuildPasses \
    +g_expectedRegexTestPasses \
    +"0) MPI_DEBUG: Will attempt to run!\n" \
    +"1) SERIAL_RELEASE: Will \*not\* attempt to run on request!\n" \
    +"0) MPI_DEBUG => passed: passed=100,notpassed=0\n" \
    +"1) SERIAL_RELEASE => Test case SERIAL_RELEASE was not run! => Does not affect push readiness!\n" \
    +g_expectedCommonOptionsSummary \
    +"=> A PUSH IS READY TO BE PERFORMED!\n" \
    +"^READY TO PUSH: Trilinos:\n" \
    ,
    \
    failRegexStrList = \
    "mailx .* trilinos-checkin-tests.*\n" \
    +"DID PUSH: Trilinos\n" \
    )


def checkin_test_configure_test(testObject, testName, optionsStr, filePassRegexStrList, \
  fileFailRegexStrList=[], modifiedFilesStr="", extraPassRegexStr="" \
  ):

  if not modifiedFilesStr:
    modifiedFilesStr = "M\tpackages/teuchos/CMakeLists.txt"

  checkin_test_run_case(
    \
    testObject,
    \
    testName,
    \
    " --allow-no-pull --configure --send-email-to= --skip-push-readiness-check" \
    +" " +optionsStr \
    ,
    \
    g_cmndinterceptsCurrentBranch \
    +"IT: eg diff --name-status origin/currentbranch; 0; '"+modifiedFilesStr+"'\n" \
    +g_cmndinterceptsConfigPasses \
    ,
    \
    True,
    \
    "Configure passed!\n" \
    +"^NOT READY TO PUSH\n" \
    +extraPassRegexStr \
    ,
    filePassRegexStrList
    ,
    fileFailRegexStrList=fileFailRegexStrList
    )


def checkin_test_configure_enables_test(testObject, testName, optionsStr, regexListStr, \
  notRegexListStr="", modifiedFilesStr="", extraPassRegexStr="" \
  ):
  checkin_test_configure_test(
     testObject,
     testName,
     "--default-builds=MPI_DEBUG "+optionsStr,
     [("MPI_DEBUG/do-configure", regexListStr)],
     [("MPI_DEBUG/do-configure", notRegexListStr)],
     modifiedFilesStr,
     extraPassRegexStr,
     )


# Used as a helper in follow-up tests  
def  g_test_ss_extra_builds_ss_do_all_pass(testObject, testName):

  testBaseDir = create_checkin_test_case_dir(testName, g_verbose)

  writeStrToFile(testBaseDir+"/COMMON.config",
    "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON\n" \
    +"-DBUILD_SHARED:BOOL=ON\n" \
    +"-DTPL_BLAS_LIBRARIES:PATH=/usr/local/libblas.a\n" \
    +"-DTPL_LAPACK_LIBRARIES:PATH=/usr/local/liblapack.a\n" \
    )

  writeStrToFile(testBaseDir+"/MPI_DEBUG_SS.config",
    "-DTPL_ENABLE_MPI:BOOL=ON\n" \
    +"-DTeuchos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON\n" \
    )

  modifiedFilesStr = ""

  checkin_test_run_case(
    \
    testObject,
    \
    testName,
    \
    "--make-options=-j3 --ctest-options=-j5" \
    +" --default-builds=MPI_DEBUG --ss-extra-builds=MPI_DEBUG_SS" \
    +" --enable-packages=Phalanx" \
    +" --do-all" \
    ,
    \
    g_cmndinterceptsCurrentBranch \
    +g_cmndinterceptsPullPasses \
    +g_cmndinterceptsSendBuildTestCaseEmail \
    +g_cmndinterceptsConfigBuildTestPasses \
    +g_cmndinterceptsSendBuildTestCaseEmail \
    +g_cmndinterceptsSendFinalEmail \
    ,
    \
    True,
    \
    "Phalanx of type SS is being excluded because it is not in the valid list of package types .PS.\n" \
    +"passed: Trilinos/MPI_DEBUG: skipped configure, build, test due to no enabled packages\n" \
    +"passed: Trilinos/MPI_DEBUG_SS: passed=100,notpassed=0\n" \
    +"0) MPI_DEBUG => Skipped configure, build, test due to no enabled packages! => Does not affect push readiness!\n" \
    +"2) MPI_DEBUG_SS => passed: passed=100,notpassed=0\n" \
    +"^READY TO PUSH\n" \
    )


#
# checkin_test unit tests
#


class test_checkin_test(unittest.TestCase):


  # A) Test basic passing use cases


  def test_help(self):
    testName = "help"
    checkin_test_run_case(
      self,
      testName,
      "--help",
      "", # No shell commands!
      True,
      "checkin-test.py \[OPTIONS\]\n" \
      "Quickstart\:\n" \
      "Detailed Documentation:\n" \
      ".*--show-defaults.*\n" \
      ,
      mustHaveCheckinTestOut=False
      )
    # Help should not write the checkin-test.out file!
    self.assertEqual(
      os.path.exists(create_checkin_test_case_dir(testName, g_verbose)+"/checkin-test.out"),
      False)


  def test_show_defaults(self):
    testName = "show_defaults"
    checkin_test_run_case(
      self,
      testName,
      "--show-defaults",
      "", # No shell commands!
      True,
      "Script: checkin-test.py\n" \
      ,
      mustHaveCheckinTestOut=False
      )
    # Help should not write the checkin-test.out file!
    self.assertEqual(
      os.path.exists(create_checkin_test_case_dir(testName, g_verbose)+"/checkin-test.out"),
      False)


  def test_do_all_push_pass(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_push_pass",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      +" --abort-gracefully-if-no-updates --do-all --push" \
      +" --execute-on-ready-to-push=\"ssh -q godel /some/dir/some_command.sh &\"",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsFinalPushPasses \
      +g_cmndinterceptsSendFinalEmail \
      +"IT: ssh -q godel /some/dir/some_command.sh &; 0; 'extra command passed'\n" \
      ,
      \
      True,
      \
      "Pulled changes from this repo!\n" \
      +"There where at least some changes pulled!\n" \
      +g_expectedRegexUpdateWithBuildCasePasses \
      +g_expectedRegexConfigPasses \
      +g_expectedRegexBuildPasses \
      +g_expectedRegexTestPasses \
      +"0) MPI_DEBUG => passed: passed=100,notpassed=0\n" \
      +"1) SERIAL_RELEASE => passed: passed=100,notpassed=0\n" \
      +g_expectedCommonOptionsSummary \
      +"=> A PUSH IS READY TO BE PERFORMED!\n" \
      +"mailx .* trilinos-checkin-tests.*\n" \
      +"^DID PUSH: Trilinos:\n" \
      +"Executing final command (ssh -q godel /some/dir/some_command.sh &) since a push is okay to be performed!\n" \
      +"Running: ssh -q godel /some/dir/some_command.sh &\n" \
      ,
      [
      (getInitialPullOutputFileName(""), "pulled changes passes\n"),
      (getModifiedFilesOutputFileName(""), "M\tpackages/teuchos/CMakeLists.txt\n"),
      (getFinalPullOutputFileName(""), "final eg pull and rebase passed\n"),
      (getFinalCommitBodyFileName(""),
         getAutomatedStatusSummaryHeaderKeyStr()+"\n"
         +"Enabled Packages: Teuchos\n" \
         +"Enabled all Forward Packages\n" \
         ),
      ("MPI_DEBUG/do-configure.base",
       "\-DTPL_ENABLE_Pthread:BOOL=OFF\n"\
       +"\-DTPL_ENABLE_BinUtils:BOOL=OFF\n"\
       +"\-DTPL_ENABLE_MPI:BOOL=ON\n" \
       +"\-DTrilinos_ENABLE_TESTS:BOOL=ON\n" \
       +"\-DCMAKE_BUILD_TYPE:STRING=RELEASE\n" \
       +"\-DTrilinos_ENABLE_DEBUG:BOOL=ON\n" \
       +"\-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON\n" \
       +"\-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON\n"),
      ("MPI_DEBUG/do-configure",
       "\./do-configure.base\n" \
       +"\-DTrilinos_ENABLE_Teuchos:BOOL=ON\n" \
       +"\-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON\n" \
       +"\-DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON\n"),
      ("SERIAL_RELEASE/do-configure.base",
       "\-DTPL_ENABLE_Pthread:BOOL=OFF\n"\
       +"\-DTPL_ENABLE_BinUtils:BOOL=OFF\n"\
       +"\-DTrilinos_ENABLE_TESTS:BOOL=ON\n" \
       +"\-DTPL_ENABLE_MPI:BOOL=OFF\n" \
       +"\-DCMAKE_BUILD_TYPE:STRING=RELEASE\n" \
       +"\-DTrilinos_ENABLE_DEBUG:BOOL=OFF\n" \
       +"\-DTrilinos_ENABLE_CHECKED_STL:BOOL=OFF\n" \
       +"\-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF\n"),
      ("SERIAL_RELEASE/do-configure",
       "\./do-configure.base\n" \
       +"\-DTrilinos_ENABLE_Teuchos:BOOL=ON\n" \
       +"\-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON\n" \
       +"\-DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON\n"),
      # ToDo: Add more files to check
      ]
      )


  # In this test, we test the behavior of the script where eg on the path
  # is not found and the default eg is used instead.  We have to test the
  # entire workflow in order to make sure that raw 'eg' is not used anywhere
  # where it matters.
  def test_do_all_no_eg_installed(self):
    eg = os.path.abspath(scriptsDir+"/../common_tools/git/eg")
    checkin_test_run_case(
      \
      self,
      \
      "do_all_no_eg_installed",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      +" --default-builds=MPI_DEBUG" \
      +" --do-all --push" \
      ,
      \
      "IT: git config --get user.email; 0; bogous@somwhere.com\n" \
      +"IT: which eg; 1; '/usr/bin/which: no eg in (path1:path2:path3)'\n" \
      +"IT: "+eg+" --version; 0; "+g_officialEgVersion+"\n" \
      +"IT: "+eg+" branch; 0; '* currentbranch'\n" \
      +"IT: "+eg+" status; 0; '(on master branch)'\n" \
      +"IT: "+eg+" pull; 0; 'initial eg pull passed'\n" \
      +"IT: "+eg+" diff --name-status origin/currentbranch; 0; 'M\tpackages/teuchos/CMakeLists.txt'\n" \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +"IT: "+eg+" pull && "+eg+" rebase --against origin/currentbranch; 0; 'final eg pull and rebase passed'\n"
      +"IT: "+eg+" cat-file -p HEAD; 0; 'This is the last commit message'\n" \
      +"IT: "+eg+" log --oneline currentbranch \^origin/currentbranch; 0; '12345 Only one commit'\n" \
      +"IT: "+eg+" log --pretty=format:'%h' currentbranch\^ \^origin/currentbranch; 0; '12345'\n"
      +"IT: "+eg+" commit --amend -F .*; 0; 'Amending the last commit passed'\n"
      +"IT: "+eg+" log --oneline currentbranch \^origin/currentbranch; 0; '54321 Only one commit'\n"
      +"IT: cat modifiedFiles.out; 0; 'M\tpackages/teuchos/CMakeLists.txt'\n"\
      +"IT: "+eg+" push; 0; 'push passes'\n" \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "Warning, the eg command is not in your path! .*no eg in .path1:path2:path3.*\n" \
      "Setting to default eg in source tree '.*/tribits/common_tools/git/eg'\n" \
      ,
      inPathEg=False, egVersion=False
      )


  def test_do_all_default_builds_mpi_debug_pass(self):
    g_test_do_all_default_builds_mpi_debug_pass(self, "do_all_default_builds_mpi_debug_pass")


  def test_local_do_all_default_builds_mpi_debug_pass(self):
    checkin_test_run_case(
      \
      self,
      \
      "local_do_all_default_builds_mpi_debug_pass",
      \
      "--make-options=-j3 --ctest-options=-j5 --default-builds=MPI_DEBUG" \
      +" --extra-pull-from=machine:/path/to/repo:master --local-do-all" \
      +" --execute-on-ready-to-push=\"ssh -q godel /some/dir/some_command.sh &\"",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      g_expectedRegexConfigPasses \
      +g_expectedRegexBuildPasses \
      +g_expectedRegexTestPasses \
      +"0) MPI_DEBUG => passed: passed=100,notpassed=0\n" \
      +"1) SERIAL_RELEASE => Test case SERIAL_RELEASE was not run! => Does not affect push readiness!\n" \
      +g_expectedCommonOptionsSummary \
      +"A current successful pull does \*not\* exist => Not ready for final push!\n" \
      +"Explanation: In order to safely push, the local working directory needs\n" \
      +"A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      +"^NOT READY TO PUSH: Trilinos:\n" \
      +"Not executing final command (ssh -q godel /some/dir/some_command.sh &) since a push is not okay to be performed!\n" \
      )


  def test_do_all_default_builds_mpi_debug_test_fail_force_push_pass(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_default_builds_mpi_debug_test_fail_force_push_pass",
      \
      "--make-options=-j3 --ctest-options=-j5 --default-builds=MPI_DEBUG" \
      " --do-all --force-push --push",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildPasses \
      +"IT: ctest -j5; 1; '80% tests passed, 20 tests failed out of 100'\n" \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsFinalPushPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,      \
      True,
      \
      g_expectedRegexUpdateWithBuildCasePasses \
      +g_expectedRegexConfigPasses \
      +g_expectedRegexBuildPasses \
      +"FAILED: ctest failed returning 1!\n" \
      +"testResultsLine = .80% tests passed, 20 tests failed out of 100.\n" \
      +"0) MPI_DEBUG => FAILED: passed=80,notpassed=20\n" \
      +"1) SERIAL_RELEASE => Test case SERIAL_RELEASE was not run! => Does not affect push readiness!\n" \
      +g_expectedCommonOptionsSummary \
      +"Test: FAILED\n" \
      +"=> A PUSH IS READY TO BE PERFORMED!\n" \
      +"\*\*\* WARNING: The acceptance criteria for doing a push has \*not\*\n" \
      +"\*\*\* been met, but a push is being forced anyway by --force-push!\n" \
      +"DID FORCED PUSH: Trilinos:\n" \
      +"REQUESTED ACTIONS: PASSED\n"
      )


  def test_do_all_default_builds_mpi_debug_then_wipe_clean_pull_pass(self):

    testName = "do_all_default_builds_mpi_debug_then_from_scratch_pull_pass"

    # Do the build/test only first (ready to push)
    g_test_do_all_default_builds_mpi_debug_pass(self, testName)

    # Do the push after the fact
    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      "--make-options=-j3 --ctest-options=-j5 --default-builds=MPI_DEBUG" \
      +" --wipe-clean --pull" \
      ,
      \
      g_cmndinterceptsCurrentBranch \
      +"FT: rm -rf MPI_DEBUG\n" \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "Running: rm -rf MPI_DEBUG\n" \
      +"0) MPI_DEBUG => No configure, build, or test for MPI_DEBUG was requested! => Not ready to push!\n" \
      +"=> A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      +"^NOT READY TO PUSH: Trilinos:\n"
      )


  def test_send_email_only_on_failure_do_all_push_pass(self):
    checkin_test_run_case(
      \
      self,
      \
      "send_email_only_on_failure_do_all_push_pass",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      +" --abort-gracefully-if-no-updates --do-all --push" \
      +" --send-email-only-on-failure" \
      ,
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsFinalPushPasses \
      ,
      \
      True,
      \
      g_expectedRegexUpdateWithBuildCasePasses \
      +g_expectedRegexConfigPasses \
      +g_expectedRegexBuildPasses \
      +g_expectedRegexTestPasses \
      +g_expectedCommonOptionsSummary \
      +"MPI_DEBUG: Skipping sending build/test case email because it passed and --send-email-only-on-failure was set!\n" \
      +"SERIAL_RELEASE: Skipping sending build/test case email because it passed and --send-email-only-on-failure was set!\n" \
      +"Skipping sending final email because it passed and --send-email-only-on-failure was set!\n" \
      )


  def test_abort_gracefully_if_no_enables(self):
    checkin_test_run_case(
      \
      self,
      \
      "abort_gracefully_if_no_enables",
      \
      " --abort-gracefully-if-no-enables --do-all --push",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsStatusPullPasses \
      +g_cmndinterceptsDiffOnlyNoChangesPasses \
      ,
      \
      True,
      \
      "SKIPPED: MPI_DEBUG configure skipped because no packages are enabled!\n" \
      "SKIPPED: MPI_DEBUG build skipped because configure did not pass!\n" \
      "SKIPPED: MPI_DEBUG tests skipped because no packages are enabled!\n" \
      "SKIPPED: SERIAL_RELEASE configure skipped because no packages are enabled!\n" \
      "SKIPPED: SERIAL_RELEASE build skipped because configure did not pass!\n" \
      "SKIPPED: SERIAL_RELEASE tests skipped because no packages are enabled!\n" \
      +"subjectLine = .passed: Trilinos/MPI_DEBUG: skipped configure, build, test due to no enabled packages.\n" \
      +"subjectLine = .passed: Trilinos/SERIAL_RELEASE: skipped configure, build, test due to no enabled packages.\n" \
      +"0) MPI_DEBUG => Skipped configure, build, test due to no enabled packages! => Does not affect push readiness!\n" \
      +"1) SERIAL_RELEASE => Skipped configure, build, test due to no enabled packages! => Does not affect push readiness!\n" \
      +"MPI_DEBUG: Skipping sending build/test case email because there were no enables and --abort-gracefully-if-no-enables was set!\n"
      +"SERIAL_RELEASE: Skipping sending build/test case email because there were no enables and --abort-gracefully-if-no-enables was set!\n"
      +"There were no successfuly attempts to configure/build/test!\n" \
      +"Skipping sending final email because there were no enables and --abort-gracefully-if-no-enables was set!\n" \
      +"ABORTED DUE TO NO ENABLES: Trilinos:\n" \
      +"REQUESTED ACTIONS: PASSED\n" \
      )


  # ToDo: Add a test case where PS has no enables but SS does!


  def test_do_all_no_append_test_results_push_pass(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_no_append_test_results_push_pass",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      +" --do-all --no-append-test-results --push",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsFinalPushNoAppendTestResultsPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      g_expectedRegexUpdateWithBuildCasePasses \
      +g_expectedRegexConfigPasses \
      +g_expectedRegexBuildPasses \
      +g_expectedRegexTestPasses \
      +"0) MPI_DEBUG => passed: passed=100,notpassed=0\n" \
      +"1) SERIAL_RELEASE => passed: passed=100,notpassed=0\n" \
      +g_expectedCommonOptionsSummary \
      +"Skipping appending test results on request (--no-append-test-results)!\n" \
      +"=> A PUSH IS READY TO BE PERFORMED!\n" \
      +"^DID PUSH: Trilinos:\n" \
      )


  def test_do_all_no_rebase_push_pass(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_no_rebase_push_pass",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      +" --do-all --no-rebase --push",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsFinalPushNoRebasePasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      g_expectedRegexUpdateWithBuildCasePasses \
      +g_expectedRegexConfigPasses \
      +g_expectedRegexBuildPasses \
      +g_expectedRegexTestPasses \
      +"0) MPI_DEBUG => passed: passed=100,notpassed=0\n" \
      +"1) SERIAL_RELEASE => passed: passed=100,notpassed=0\n" \
      +g_expectedCommonOptionsSummary \
      +"Skipping the final rebase on request! (see --no-rebase option)\n" \
      +"=> A PUSH IS READY TO BE PERFORMED!\n" \
      +"^DID PUSH: Trilinos:\n" \
      )


  def test_extra_repo_1_explicit_enable_configure_pass(self):

    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"

    testName = "extra_repo_1_explicit_enable_configure_pass"

    testBaseDir = create_checkin_test_case_dir(testName, g_verbose)

    writeStrToFile(testBaseDir+"/MPI_DEBUG_SS.config",
      "-DTPL_ENABLE_MPI:BOOL=ON\n" \
      +"-DTeuchos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON\n" \
      )

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      "--extra-repos=preCopyrightTrilinos --allow-no-pull --without-default-builds" \
      " --extra-builds=MPI_DEBUG_SS --enable-packages=Stalix --configure", \
      \
      "IT: .*cmake .+ -P .+/TribitsDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
      +g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsDiffOnlyPassesPreCopyrightTrilinos \
      +g_cmndinterceptsConfigPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "-extra-repos=.preCopyrightTrilinos.\n" \
      +"Pulling in packages from extra repos: preCopyrightTrilinos ...\n" \
      +"projectDepsXmlFileOverride="+projectDepsXmlFileOverride+"\n" \
      +"Enabling only the explicitly specified packages .Stalix. ...\n" \
      +"Trilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos\n" \
      +"Enabled Packages: Stalix\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_implicit_enable_configure_pass(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_implicit_enable_configure_pass",
      \
      "--extra-repos=preCopyrightTrilinos --allow-no-pull --default-builds=MPI_DEBUG --configure", \
      \
      "IT: .*cmake .+ -P .+/TribitsDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
      +g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsDiffOnlyPassesPreCopyrightTrilinos \
      +g_cmndinterceptsConfigPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "-extra-repos=.preCopyrightTrilinos.\n" \
      +"Pulling in packages from extra repos: preCopyrightTrilinos ...\n" \
      +"projectDepsXmlFileOverride="+projectDepsXmlFileOverride+"\n" \
      +"Modified file: .preCopyrightTrilinos/teko/CMakeLists.txt.\n" \
      +"  => Enabling .Teko.!\n" \
      +"Teko of type SS is being excluded because it is not in the valid list of package types .PS.\n" \
      +"Trilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos\n" \
      +"Enabled Packages: Teuchos, Teko\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_do_all_push_pass(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_do_all_push_pass",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --default-builds=MPI_DEBUG --do-all --push", \
      \
      g_cmndinterceptsExtraRepo1DoAllUpToPush \
      +g_cmndinterceptsCatModifiedFilesPasses \
      +g_cmndinterceptsPushOnlyPasses \
      +g_cmndinterceptsCatModifiedFilesPreCoprightTrilinosPasses \
      +g_cmndinterceptsPushOnlyPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "-extra-repos=.preCopyrightTrilinos.\n" \
      +"Pulling in packages from extra repos: preCopyrightTrilinos ...\n" \
      +"Enabling .Teko..\n" \
      +"Teko of type SS is being excluded because it is not in the valid list of package types .PS.\n" \
      +"pullInitial.preCopyrightTrilinos.out\n" \
      +"Update passed!\n"\
      +"All of the tests ran passed!\n" \
      +"pullFinal.preCopyrightTrilinos.out\n" \
      +"Final update passed!\n" \
      +"commitFinalBody.preCopyrightTrilinos.out\n" \
      +"commitFinal.preCopyrightTrilinos.out\n" \
      +"push.preCopyrightTrilinos.out\n" \
      +"Push passed!\n" \
      +"Enabled Packages: Teuchos, Teko\n" \
      +"DID PUSH: Trilinos:\n" \
      +"REQUESTED ACTIONS: PASSED\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_pull_extra_pull_pass(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_pull_extra_pull_pass",
      \
      "--extra-repos=preCopyrightTrilinos --pull --extra-pull-from=somemachine:someotherbranch", \
      \
      "IT: .*cmake .+ -P .+/TribitsDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
      +g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +"IT: eg pull somemachine someotherbranch; 0; 'eg extra pull passed'\n"
      +"IT: eg pull somemachine someotherbranch; 0; 'eg extra pull passed'\n"
      +g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsDiffOnlyPassesPreCopyrightTrilinos \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "pullInitial.out\n" \
      "pullInitial.preCopyrightTrilinos.out\n" \
      "pullInitialExtra.out\n" \
      "pullInitialExtra.preCopyrightTrilinos.out\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_trilinos_changes_do_all_push_pass(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_trilinos_changes_do_all_push_pass",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --default-builds=MPI_DEBUG --do-all --push", \
      \
      g_cmndinterceptsExtraRepo1DoAllUpToPush \
      +g_cmndinterceptsCatModifiedFilesPasses \
      +g_cmndinterceptsPushOnlyPasses \
      +g_cmndinterceptsCatModifiedFilesPreCoprightTrilinosNoChanges \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "Skipping push to .preCopyrightTrilinos. because there are no changes!\n" \
      +"Push passed!\n" \
      +"DID PUSH: Trilinos:\n" \
      +"REQUESTED ACTIONS: PASSED\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_extra_repo_changes_do_all_push_pass(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_extra_repo_changes_do_all_push_pass",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --default-builds=MPI_DEBUG --do-all --push", \
      \
      g_cmndinterceptsExtraRepo1DoAllUpToPush \
      +g_cmndinterceptsCatModifiedFilesNoChanges \
      +g_cmndinterceptsCatModifiedFilesPreCoprightTrilinosPasses \
      +g_cmndinterceptsPushOnlyPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "Skipping push to .. because there are no changes!\n" \
      +"Push passed!\n" \
      +"DID PUSH: Trilinos:\n" \
      +"REQUESTED ACTIONS: PASSED\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_abort_gracefully_if_no_updates_no_updates_passes(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_abort_gracefully_if_no_updates_no_updates_passes",
      \
      "--extra-repos=preCopyrightTrilinos --abort-gracefully-if-no-updates --do-all --pull", \
      \
      g_cmndinterceptsExtraRepo1ThroughStatusPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsDiffOnlyNoChangesPasses \
      +g_cmndinterceptsDiffOnlyPassesPreCopyrightTrilinos \
      ,
      \
      True,
      \
      "Pulling in packages from extra repos: preCopyrightTrilinos ...\n" \
      +"Did not pull any changes from this repo!\n" \
      +"No changes were pulled!\n" \
      +"Not perfoming any build cases because pull did not give any changes" \
        " and --abort-gracefully-if-no-updates!\n" \
      +"Skipping sending final email because there were no updates" \
          " and --abort-gracefully-if-no-updates was set!\n" \
      +"ABORTED DUE TO NO UPDATES\n" \
      +"REQUESTED ACTIONS: PASSED\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_extra_pull_abort_gracefully_if_no_updates_no_updates_passes(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_extra_pull_abort_gracefully_if_no_updates_no_updates_passes",
      \
      "--extra-repos=preCopyrightTrilinos --abort-gracefully-if-no-updates" \
      +" --extra-pull-from=machine:master --do-all --pull" \
      ,
      \
      g_cmndinterceptsExtraRepo1ThroughStatusPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsDiffOnlyNoChangesPasses \
      +g_cmndinterceptsDiffOnlyPassesPreCopyrightTrilinos \
      ,
      \
      True,
      \
      "Pulling in packages from extra repos: preCopyrightTrilinos ...\n" \
      +"Did not pull any changes from this repo!\n" \
      +"No changes were pulled!\n" \
      +"Not perfoming any build cases because pull did not give any changes" \
        " and --abort-gracefully-if-no-updates!\n" \
      +"Skipping sending final email because there were no updates" \
          " and --abort-gracefully-if-no-updates was set!\n" \
      +"ABORTED DUE TO NO UPDATES\n" \
      +"REQUESTED ACTIONS: PASSED\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_extra_pull_abort_gracefully_if_no_updates_main_repo_update(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_extra_pull_abort_gracefully_if_no_updates_main_repo_update",
      \
      "--extra-repos=preCopyrightTrilinos --abort-gracefully-if-no-updates" \
      +" --extra-pull-from=machine:master --pull" \
      ,
      \
      g_cmndinterceptsExtraRepo1ThroughStatusPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsDiffOnlyNoChangesPasses \
      +g_cmndinterceptsDiffOnlyPassesPreCopyrightTrilinos \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "Pulled changes from this repo!\n" \
      +"Did not pull any changes from this repo!\n" \
      +"There where at least some changes pulled!\n" \
      +"Update passed!\n" \
      +"NOT READY TO PUSH\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_extra_pull_abort_gracefully_if_no_updates_extra_repo_update(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_extra_pull_abort_gracefully_if_no_updates_extra_repo_update",
      \
      "--extra-repos=preCopyrightTrilinos --abort-gracefully-if-no-updates" \
      +" --extra-pull-from=machine:master --pull" \
      ,
      \
      g_cmndinterceptsExtraRepo1ThroughStatusPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsDiffOnlyNoChangesPasses \
      +g_cmndinterceptsDiffOnlyPassesPreCopyrightTrilinos \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "Pulled changes from this repo!\n" \
      +"Did not pull any changes from this repo!\n" \
      +"There where at least some changes pulled!\n" \
      +"Update passed!\n" \
      +"NOT READY TO PUSH\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_extra_pull_abort_gracefully_if_no_updates_main_repo_extra_update(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_extra_pull_abort_gracefully_if_no_updates_main_repo_extra_update",
      \
      "--extra-repos=preCopyrightTrilinos --abort-gracefully-if-no-updates" \
      +" --extra-pull-from=machine:master --pull" \
      ,
      \
      g_cmndinterceptsExtraRepo1ThroughStatusPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsDiffOnlyNoChangesPasses \
      +g_cmndinterceptsDiffOnlyPassesPreCopyrightTrilinos \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "Pulled changes from this repo!\n" \
      +"Did not pull any changes from this repo!\n" \
      +"There where at least some changes pulled!\n" \
      +"Update passed!\n" \
      +"NOT READY TO PUSH\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_extra_pull_abort_gracefully_if_no_updates_extra_repo_extra_update(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_extra_pull_abort_gracefully_if_no_updates_extra_repo_extra_update",
      \
      "--extra-repos=preCopyrightTrilinos --abort-gracefully-if-no-updates" \
      +" --extra-pull-from=machine:master --pull" \
      ,
      \
      g_cmndinterceptsExtraRepo1ThroughStatusPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsPullOnlyNoUpdatesPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsDiffOnlyNoChangesPasses \
      +g_cmndinterceptsDiffOnlyPassesPreCopyrightTrilinos \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "Pulled changes from this repo!\n" \
      +"Did not pull any changes from this repo!\n" \
      +"There where at least some changes pulled!\n" \
      +"Update passed!\n" \
      +"NOT READY TO PUSH\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_file_2_continuous_pull(self):

    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"

    testName = "test_extra_repo_file_2_continuous_pull"

    testBaseDir = create_checkin_test_case_dir(testName, g_verbose)

    writeStrToFile(testBaseDir+"/MPI_DEBUG_SS.config",
      "-DTPL_ENABLE_MPI:BOOL=ON\n" \
      +"-DTeuchos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON\n" \
      )

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      " --extra-repos-file="+scriptsDir+"/UnitTests/ExtraReposListExisting_2.cmake" \
      " --extra-repos-type=Continuous" \
      " --extra-builds=MPI_DEBUG_SS --enable-packages=Stalix --pull", \
      \
      "IT: .*cmake .+ -P .+/TribitsDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
      +g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsDiffOnlyPassesPreCopyrightTrilinos \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "-extra-repos-file=.*ExtraReposListExisting_2.\n" \
      +"-extra-repos-type=.Continuous.\n" \
      +"Pulling in packages from extra repos: preCopyrightTrilinos ...\n" \
      +"projectDepsXmlFileOverride="+projectDepsXmlFileOverride+"\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_file_2_nightly_pull(self):

    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.extraTrilinosRepo.gold.xml"

    testName = "test_extra_repo_file_2_nightly_pull"

    testBaseDir = create_checkin_test_case_dir(testName, g_verbose)

    writeStrToFile(testBaseDir+"/MPI_DEBUG_SS.config",
      "-DTPL_ENABLE_MPI:BOOL=ON\n" \
      +"-DTeuchos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON\n" \
      )

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      " --extra-repos-file="+scriptsDir+"/UnitTests/ExtraReposListExisting_2.cmake" \
      " --extra-repos-type=Nightly" \
      " --extra-builds=MPI_DEBUG_SS --pull", \
      \
      "IT: .*cmake .+ -P .+/TribitsDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
      +g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsDiffOnlyPassesPreCopyrightTrilinos \
      +g_cmndinterceptsDiffOnlyPassesExtraTrilinosRepo \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "-extra-repos-file=.*ExtraReposListExisting_2.\n" \
      +"-extra-repos-type=.Nightly.\n" \
      +"Pulling in packages from extra repos: preCopyrightTrilinos,extraTrilinosRepo ...\n" \
      +"cmake .* -DTrilinos_EXTRA_REPOSITORIES=.preCopyrightTrilinos.extraTrilinosRepo. -P .*TribitsDumpDepsXmlScript.cmake\n" \
      +"projectDepsXmlFileOverride="+projectDepsXmlFileOverride+"\n" \
      +"3.a.1) Git Repo: .preCopyrightTrilinos.\n" \
      +"3.a.2) Git Repo: .extraTrilinosRepo.\n" \
      +"Running in working directory: .*MockTrilinos/preCopyrightTrilinos ...\n" \
      +"Running in working directory: .*MockTrilinos/extraTrilinosRepo ...\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_file_3_continuous_pull_configure(self):

    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"

    testName = "test_extra_repo_file_3_continuous_pull_configure"

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      " --extra-repos-file="+scriptsDir+"/UnitTests/ExtraReposListExisting_3.cmake" \
      " --extra-repos-type=Continuous" \
      " --default-builds=MPI_DEBUG --pull --configure", \
      \
      "IT: .*cmake .+ -P .+/TribitsDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
      +g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +"IT: eg diff --name-status origin/currentbranch; 0; ''\n" \
      +"IT: eg diff --name-status origin/currentbranch; 0; ''\n" \
      +"IT: eg diff --name-status origin/currentbranch; 0; 'M\tExtraTeuchosStuff.hpp'\n" \
      +g_cmndinterceptsConfigPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "-extra-repos-file=.*ExtraReposListExisting_3.\n" \
      +"-extra-repos-type=.Continuous.\n" \
      +"Pulling in packages from extra repos: preCopyrightTrilinos ...\n" \
      +"projectDepsXmlFileOverride="+projectDepsXmlFileOverride+"\n" \
      +"cmake .* -DTrilinos_EXTRA_REPOSITORIES=.preCopyrightTrilinos. -P .*TribitsDumpDepsXmlScript.cmake\n" \
      +"Modified file: .packages/teuchos/extrastuff/ExtraTeuchosStuff.hpp.\n" \
      +"=> Enabling .Teuchos.!\n" \
      +"Full package enable list: .Teuchos.\n" \
      ,
      [
      ("MPI_DEBUG/do-configure",
       "\-DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos\n"),
      ] \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_file_3_continuous_do_all_push(self):

    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.gold.xml"

    testName = "test_extra_repo_file_3_continuous_do_all_push"

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      " --extra-repos-file="+scriptsDir+"/UnitTests/ExtraReposListExisting_3.cmake" \
      " --extra-repos-type=Continuous" \
      " --extra-repos=ExtraTeuchosRepo" \
      " --make-options=-j3 --ctest-options=-j5" \
      " --default-builds=MPI_DEBUG --do-all --push", \
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +"IT: eg diff --name-status origin/currentbranch; 0; ''\n" \
      +"IT: eg diff --name-status origin/currentbranch; 0; 'M\tExtraTeuchosStuff.hpp'\n" \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsFinalPullRebasePasses \
      +g_cmndinterceptsFinalPullRebasePasses \
      +g_cmnginterceptsEgLogCmnds \
      +g_cmndinterceptsAmendCommitPasses \
      +g_cmnginterceptsEgLogCmnds \
      +g_cmndinterceptsAmendCommitPasses \
      +g_cmndinterceptsLogCommitsPasses \
      +g_cmndinterceptsLogCommitsPasses \
      +"IT: cat modifiedFiles.out; 0; ''\n"\
      +"IT: cat modifiedFiles.ExtraTeuchosRepo.out; 0; 'M\tExtraTeuchosStuff.hpp'\n"\
      +"IT: eg push; 0; 'push passes'\n" \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "-extra-repos-file=.*ExtraReposListExisting_3.\n" \
      +"-extra-repos-type=.Continuous.\n" \
      +"Modified file: .packages/teuchos/extrastuff/ExtraTeuchosStuff.hpp.\n" \
      +"=> Enabling .Teuchos.!\n" \
      +"Full package enable list: .Teuchos.\n" \
      +"Skipping push to .. because there are no changes!\n" \
      +"push.ExtraTeuchosRepo.out\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_file_project_nightly_nothing_fail(self):

    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"

    testName = "test_extra_repo_file_project_nightly_nothing_fail"

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      " --extra-repos-file=project" \
      " --extra-repos-type=Nightly" , \
      \
      "" \
      ,
      \
      False,
      \
      "ERROR! Skipping missing extra repo .Dakota. since\n" \
      "MockTrilinos/packages/TriKota/Dakota\n" \
      "Error, the command ..*cmake .*TribitsGetExtraReposForCheckinTest.cmake\n" \
      ,
      mustHaveCheckinTestOut=False
      )


  def test_extra_repo_file_project_continuous_extra_repos_pull(self):

    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"

    testName = "test_extra_repo_file_project_continuous_extra_repos_pull"

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      " --extra-repos-file=project" \
      " --extra-repos-type=Continuous" \
      " --extra-repos=preCopyrightTrilinos" \
      " --pull", \
      \
      "IT: .*cmake .+ -P .+/TribitsDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
      +g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsDiffOnlyPassesPreCopyrightTrilinos \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "-extra-repos-file=.project.\n" \
      +"-extra-repos-type=.Continuous.\n" \
      +"Pulling in packages from extra repos: preCopyrightTrilinos ...\n" \
      +"projectDepsXmlFileOverride="+projectDepsXmlFileOverride+"\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_file_missing_assert_fail(self):

    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"

    testName = "test_extra_repo_file_missing_assert_fail"

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      " --extra-repos-file="+scriptsDir+"/UnitTests/ExtraReposListExisting1Missing1.cmake" \
      " --extra-repos-type=Continuous" , \
      \
      "" \
      ,
      \
      False,
      \
      "ERROR! Skipping missing extra repo .MissingRepo. since\n" \
      "Error, the command ..*cmake .*TribitsGetExtraReposForCheckinTest.cmake\n" \
      ,
      mustHaveCheckinTestOut=False
      )



  def test_extra_repo_file_missing_ignore_pull(self):

    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"

    testName = "test_extra_repo_file_missing_ignore_pull"

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      " --extra-repos-file="+scriptsDir+"/UnitTests/ExtraReposListExisting1Missing1.cmake" \
      " --extra-repos-type=Continuous --ignore-missing-extra-repos --pull" , \
      \
      "IT: .*cmake .+ -P .+/TribitsDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
      +g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsDiffOnlyPassesPreCopyrightTrilinos \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "WARNING!  Ignoring missing extra repo .MissingRepo. as requested since\n" \
      "Pulling in packages from extra repos: preCopyrightTrilinos ...\n" \
      ,
      mustHaveCheckinTestOut=False
      )


#  def test_extra_repo_file_default_pull_configure(self):
#    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
#    checkin_test_run_case(
#      \
#      self,
#      \
#      "extra_repo_file_default_pull_configure",
#      \
#      "--extra-repos-file=default --pull --configure", \
#      \
#      "IT: .*cmake .+ -P .+/TribitsDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
#      +g_cmndinterceptsCurrentBranch \
#      +g_cmndinterceptsDiffOnlyPasses \
#      ,
#      \
#      True,
#      \
#      "-extra-repos-file=.default.\n"+ \
#      "Dummy" \
#      ,
#      \
#      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
#      )


  def test_abort_gracefully_if_no_updates_status_fails(self):
    checkin_test_run_case(
      \
      self,
      \
      "abort_gracefully_if_no_updates_status_fails",
      \
      "--abort-gracefully-if-no-updates --do-all --pull" \
      ,
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsStatusChangedButNotUpdatedPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      "ERROR: There are changed unstaged uncommitted files => cannot continue!\n" \
      +"No changes were pulled!\n" \
      +"Skipping getting list of modified files because pull failed!\n" \
      +"Not running any build/test cases because the update (pull) failed!\n" \
      +"  => A PUSH IS .NOT. READY TO BE PERFORMED!\n" \
      +"INITIAL PULL FAILED\n" \
      +"REQUESTED ACTIONS: FAILED\n" \
      )

    
  # B) Test package enable/disable logic


  # NOTE: The setting of built-in cmake cache variables in do-configure[.base]
  # files is tested in the unit test test_do_all_commit_push_pass(...)


  def test_read_config_files_mpi_debug(self):
    
    testName = "read_config_files"

    testBaseDir = create_checkin_test_case_dir(testName, g_verbose)

    writeStrToFile(testBaseDir+"/COMMON.config",
      "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON\n" \
      +"-DBUILD_SHARED:BOOL=ON\n" \
      +"-DTPL_BLAS_LIBRARIES:PATH=/usr/local/libblas.a\n" \
      +"-DTPL_LAPACK_LIBRARIES:PATH=/usr/local/liblapack.a\n" \
      )

    writeStrToFile(testBaseDir+"/MPI_DEBUG.config",
      "-DMPI_BASE_DIR:PATH=/usr/lib64/openmpi/1.2.7-gcc\n" \
      "-DMPI_CXX_COMPILER:PATHNAME=/usr/lib64/openmpi/1.2.7-gcc/mpicxx\n" \
      "-DMPI_C_COMPILER:PATHNAME=/usr/lib64/openmpi/1.2.7-gcc/mpicc\n" \
      "-DMPI_Fortran_COMPILER:PATHNAME=/usr/lib64/openmpi/1.2.7-gcc/mpif77\n" \
      )

    checkin_test_configure_test(
      \
      self,
      \
      testName,
      \
      "--default-builds=MPI_DEBUG",
      \
      [
      ("MPI_DEBUG/do-configure.base",
       "\-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON\n" \
       +"\-DBUILD_SHARED:BOOL=ON\n" \
       +"\-DTPL_BLAS_LIBRARIES:PATH=/usr/local/libblas.a\n" \
       +"\-DTPL_LAPACK_LIBRARIES:PATH=/usr/local/liblapack.a\n"
       +"\-DMPI_BASE_DIR:PATH=/usr/lib64/openmpi/1.2.7-gcc\n" \
       +"\-DMPI_CXX_COMPILER:PATHNAME=/usr/lib64/openmpi/1.2.7-gcc/mpicxx\n" \
       +"\-DMPI_C_COMPILER:PATHNAME=/usr/lib64/openmpi/1.2.7-gcc/mpicc\n" \
       +"\-DMPI_Fortran_COMPILER:PATHNAME=/usr/lib64/openmpi/1.2.7-gcc/mpif77\n" \
       ),
      ]
      )


  def test_auto_enable(self):
    checkin_test_configure_enables_test(
      \
      self,
      \
      "auto_enable",
      \
      "", # Allow auto-enable of Teuchos!
      \
      "\-DTrilinos_ENABLE_Teuchos:BOOL=ON\n" \
      )


  def test_enable_packages(self):
    checkin_test_configure_enables_test(
      \
      self,
      \
      "enable_packages",
      \
      "--enable-packages=TrilinosFramework,RTOp,Thyra",
      \
      "\-DTrilinos_ENABLE_TrilinosFramework:BOOL=ON\n" \
      +"\-DTrilinos_ENABLE_RTOp:BOOL=ON\n" \
      +"\-DTrilinos_ENABLE_Thyra:BOOL=ON\n" \
      ,
      \
      "\-DTrilinos_ENABLE_Teuchos:BOOL=ON\n" \
      )
    # Above, the --enable-packages option turns off the check of the modified
    # files and set the enables manually.


  def test_disable_packages(self):
    checkin_test_configure_enables_test(
      \
      self,
      \
      "disable_packages",
      \
      "--disable-packages=Tpetra,Thyra",
      \
      "\-DTrilinos_ENABLE_Teuchos:BOOL=ON\n" \
      +"\-DTrilinos_ENABLE_Tpetra:BOOL=OFF\n" \
      +"\-DTrilinos_ENABLE_Thyra:BOOL=OFF\n" \
      )
    # Above: --disable-packages does not turn off auto-enable and therefore
    # Teuchos is picked up.


  def test_enable_disable_packages(self):
    checkin_test_configure_enables_test(
      \
      self,
      \
      "enable_disable_packages",
      \
      "--enable-packages=TrilinosFramework,RTOp,Thyra,Tpetra" \
      +" --disable-packages=Tpetra,Stratimikos",
      \
      "\-DTrilinos_ENABLE_TrilinosFramework:BOOL=ON\n" \
      +"\-DTrilinos_ENABLE_RTOp:BOOL=ON\n" \
      +"\-DTrilinos_ENABLE_Thyra:BOOL=ON\n" \
      +"\-DTrilinos_ENABLE_Tpetra:BOOL=OFF\n" \
      +"\-DTrilinos_ENABLE_Stratimikos:BOOL=OFF\n" \
      ,
      \
      "\-DTrilinos_ENABLE_Teuchos:BOOL=ON\n" \
      +"\-DTrilinos_ENABLE_Tpetra:BOOL=ON\n" \
      )
    # Above, Teuchos should not be enabled because --enable-packages should
    # result in the modified file in Teuchos to be ignored.  The enable for
    # Tpetra should not be on because it should be removed from the enable
    # list.


  def test_no_enable_fwd_packages(self):
    checkin_test_configure_enables_test(
      self,
      "no_enable_fwd_packages",
      "--no-enable-fwd-packages",
      "\-DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=OFF\n" \
      )


  def test_enable_all_packages_auto_implicit(self):
    checkin_test_configure_enables_test(
      self,
      "enable_all_packages_auto",
      "", # --enable-all-packages=auto
      "\-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON\n" \
      +"\-DTrilinos_ENABLE_TrilinosFramework:BOOL=ON\n",
      modifiedFilesStr="M\tcmake/utils/AppendSet.cmake",
      )


  def test_enable_all_packages_auto(self):
    checkin_test_configure_enables_test(
      self,
      "enable_all_packages_auto",
      "--enable-all-packages=auto",
      "\-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON\n",
      modifiedFilesStr="M\tCMakeLists.txt", # Will not trigger TrilinosFramework!
      extraPassRegexStr="Modifed file: .CMakeLists.txt.\n"\
      +"Enabling all Trilinos packages!\n",
      )


  def test_enable_all_packages_on(self):
    checkin_test_configure_enables_test(
      self,
      "enable_all_packages_on",
      "--enable-all-packages=on",
      "\-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON\n",
      modifiedFilesStr = "M\tdummy.txt", # Will not trigger any enables!
      extraPassRegexStr="Enabling all packages on request\n",
      )


  def test_enable_all_packages_off(self):
    checkin_test_configure_enables_test(
      self,
      "enable_all_packages_auto",
      "--enable-all-packages=off",
      "\-DTrilinos_ENABLE_TrilinosFramework:BOOL=ON\n",
      notRegexListStr="\-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON\n",
      modifiedFilesStr="M\tcmake/utils/AppendSet.cmake",
      )


  # C) Test partial actions short of running tests


  def test_default_builds_mpi_debug_pull_only(self):
    checkin_test_run_case(
      self,
      \
      "default_builds_mpi_debug_pull_only",
      \
      "--default-builds=MPI_DEBUG --pull",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      g_expectedRegexUpdatePasses \
      +"Not performing any build cases because no --configure, --build or --test was specified!\n" \
      +"0) MPI_DEBUG => No configure, build, or test for MPI_DEBUG was requested! => Not ready to push!\n" \
      +"A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      +"^NOT READY TO PUSH: Trilinos:\n"
      )


  def test_default_builds_mpi_debug_pull_skip_push_readiness_check(self):
    checkin_test_run_case(
      self,
      \
      "default_builds_mpi_debug_pull_skip_push_readiness_check",
      \
      "--default-builds=MPI_DEBUG --pull --skip-push-readiness-check",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      ,
      \
      True,
      \
      g_expectedRegexUpdatePasses \
      +"Skipping push readiness check on request!\n" \
      +"Not performing push or sending out push readiness status on request!\n" \
      "^NOT READY TO PUSH$\n" \
      +"REQUESTED ACTIONS: PASSED\n"
      )


  def test_default_builds_mpi_debug_pull_extra_pull_only(self):
    checkin_test_run_case(
      self,
      \
      "default_builds_mpi_debug_pull_extra_pull_only",
      \
      "--pull --extra-pull-from=machine:/repo/dir/repo:master",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsStatusPullPasses \
      +"IT: eg pull machine:/repo/dir/repo master; 0; 'eg extra pull passed'\n"
      +g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      g_expectedRegexUpdatePasses \
      +"Pulling in updates from .machine:\/repo\/dir\/repo master.\n" \
      +"eg pull machine:\/repo\/dir\/repo master\n" \
      +"Not performing any build cases because no --configure, --build or --test was specified!\n" \
      +"A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      +"^NOT READY TO PUSH: Trilinos:\n"
      )


  def test_default_builds_mpi_debug_extra_pull_only(self):
    checkin_test_run_case(
      self,
      \
      "default_builds_mpi_debug_extra_pull_only",
      \
      "--extra-pull-from=machine:master",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      "Skipping all updates on request!\n" \
      +"Not performing any build cases because no --configure, --build or --test was specified!\n" \
      +"A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      +"^INITIAL PULL FAILED: Trilinos:\n"
      )


  def test_default_builds_mpi_debug_configure_only(self):
    checkin_test_run_case(
      self,
      \
      "default_builds_mpi_debug_configure_only",
      \
      "--default-builds=MPI_DEBUG --pull --configure",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      g_expectedRegexUpdateWithBuildCasePasses+ \
      "Configure passed!\n" \
      "touch configure.success\n" \
      "Skipping the build on request!\n" \
      "Skipping the tests on request!\n" \
      "0) MPI_DEBUG => passed: configure-only passed => Not ready to push!\n" \
      "A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      "NOT READY TO PUSH: Trilinos:\n"
      )


  def test_default_builds_mpi_debug_build_only(self):
    checkin_test_run_case(
      self,
      \
      "default_builds_mpi_debug_build_only",
      \
      "--make-options=-j3 --default-builds=MPI_DEBUG --pull --configure --build",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      g_expectedRegexUpdateWithBuildCasePasses+ \
      "Configure passed!\n" \
      "touch configure.success\n" \
      "Build passed!\n" \
      "Skipping the tests on request!\n" \
      "0) MPI_DEBUG => passed: build-only passed => Not ready to push!\n" \
      "A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      "NOT READY TO PUSH: Trilinos:\n"
      )


  # D) Test --extra-builds and --ss-extra-builds


  def test_extra_builds_read_config_file(self):
    
    testName = "extra_builds_read_config_file"

    testBaseDir = create_checkin_test_case_dir(testName, g_verbose)

    writeStrToFile(testBaseDir+"/COMMON.config",
      "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON\n" \
      +"-DBUILD_SHARED:BOOL=ON\n" \
      +"-DTPL_BLAS_LIBRARIES:PATH=/usr/local/libblas.a\n" \
      +"-DTPL_LAPACK_LIBRARIES:PATH=/usr/local/liblapack.a\n" \
      )

    writeStrToFile(testBaseDir+"/SERIAL_DEBUG_BOOST_TRACING.config",
      "-DTPL_ENABLE_BOOST:BOOL=ON\n" \
      +"-DTeuchos_ENABLE_RCP_NODE_TRACING:BOOL=ON\n" \
      )

    checkin_test_configure_test(
      \
      self,
      \
      testName,
      \
      "--without-default-builds --extra-builds=SERIAL_DEBUG_BOOST_TRACING",
      \
      [
      ("SERIAL_DEBUG_BOOST_TRACING/do-configure.base",
       "\-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON\n" \
       +"\-DBUILD_SHARED:BOOL=ON\n" \
       +"\-DTPL_BLAS_LIBRARIES:PATH=/usr/local/libblas.a\n" \
       +"\-DTPL_LAPACK_LIBRARIES:PATH=/usr/local/liblapack.a\n"
       +"\-DTPL_ENABLE_BOOST:BOOL=ON\n" \
       +"\-DTeuchos_ENABLE_RCP_NODE_TRACING:BOOL=ON\n" \
       ),
      ]
      )


  def test_extra_builds_missing_config_file_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "extra_builds_missing_config_file_fail",
      \
      "--extra-builds=SERIAL_DEBUG_BOOST_TRACING",
      \
      "", # No shell commands!
      \
      False,
      \
      "Error, the extra build configuration file SERIAL_DEBUG_BOOST_TRACING.config" \
      +" does not exit!\n" \
      )


  def test_ss_extra_builds_ps_only_pass(self):
    
    testName = "ss_extra_builds_ps_only_pass"

    testBaseDir = create_checkin_test_case_dir(testName, g_verbose)

    writeStrToFile(testBaseDir+"/COMMON.config",
      "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON\n" \
      +"-DBUILD_SHARED:BOOL=ON\n" \
      +"-DTPL_BLAS_LIBRARIES:PATH=/usr/local/libblas.a\n" \
      +"-DTPL_LAPACK_LIBRARIES:PATH=/usr/local/liblapack.a\n" \
      )

    writeStrToFile(testBaseDir+"/MPI_DEBUG_SS.config",
      "-DTPL_ENABLE_MPI:BOOL=ON\n" \
      +"-DTeuchos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON\n" \
      )

    modifiedFilesStr = "M\tpackages/teuchos/CMakeLists.txt"

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      "--make-options=-j3 --ctest-options=-j5" \
      +" --default-builds=MPI_DEBUG --do-all --push " \
      +" --ss-extra-builds=MPI_DEBUG_SS" \
      ,
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsFinalPushPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
      +"passed: Trilinos/MPI_DEBUG_SS: passed=100,notpassed=0\n" \
      +"0) MPI_DEBUG => passed: passed=100,notpassed=0\n" \
      +"2) MPI_DEBUG_SS => passed: passed=100,notpassed=0\n" \
      +"^DID PUSH\n" \
      )


  def test_ss_extra_builds_skip_case_no_email_ex_only_pass(self):
    
    testName = "ss_extra_builds_skip_case_no_email_ex_only_pass"

    testBaseDir = create_checkin_test_case_dir(testName, g_verbose)

    writeStrToFile(testBaseDir+"/MPI_DEBUG_SS.config",
      "-DTPL_ENABLE_MPI:BOOL=ON\n" \
      )

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      "--make-options=-j3 --ctest-options=-j5" \
      +" --default-builds=MPI_DEBUG" \
      +" --skip-case-no-email --do-all --push " \
      +" --extra-builds=MPI_DEBUG_SS" \
      ,
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +"IT: eg diff --name-status origin/currentbranch; 0; 'M\tpackages/stokhos/CMakeLists.txt'\n" \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsFinalPushPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "Skipping sending final status email for MPI_DEBUG because it had no packages enabled and --skip-case-no-email was set!\n" \
      +"^DID PUSH\n" \
      )


  def test_ss_extra_builds_ss_only_pass(self):
    
    testName = "ss_extra_builds_ss_only_pass"

    testBaseDir = create_checkin_test_case_dir(testName, g_verbose)

    writeStrToFile(testBaseDir+"/COMMON.config",
      "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON\n" \
      +"-DBUILD_SHARED:BOOL=ON\n" \
      +"-DTPL_BLAS_LIBRARIES:PATH=/usr/local/libblas.a\n" \
      +"-DTPL_LAPACK_LIBRARIES:PATH=/usr/local/liblapack.a\n" \
      )

    writeStrToFile(testBaseDir+"/MPI_DEBUG_SS.config",
      "-DTPL_ENABLE_MPI:BOOL=ON\n" \
      +"-DTeuchos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON\n" \
      )

    modifiedFilesStr = ""

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      "--make-options=-j3 --ctest-options=-j5" \
      +" --default-builds=MPI_DEBUG --do-all --push " \
      +" --ss-extra-builds=MPI_DEBUG_SS --enable-packages=Phalanx" \
      ,
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsFinalPushPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "Phalanx of type SS is being excluded because it is not in the valid list of package types .PS.\n" \
      +"passed: Trilinos/MPI_DEBUG: skipped configure, build, test due to no enabled packages\n" \
      +"passed: Trilinos/MPI_DEBUG_SS: passed=100,notpassed=0\n" \
      +"0) MPI_DEBUG => Skipped configure, build, test due to no enabled packages! => Does not affect push readiness!\n" \
      +"2) MPI_DEBUG_SS => passed: passed=100,notpassed=0\n" \
      +"^DID PUSH\n" \
      )


  def test_ss_extra_builds_ps_ss_pass(self):
    
    testName = "ss_extra_builds_ps_ss_pass"

    testBaseDir = create_checkin_test_case_dir(testName, g_verbose)

    writeStrToFile(testBaseDir+"/COMMON.config",
      "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON\n" \
      +"-DBUILD_SHARED:BOOL=ON\n" \
      +"-DTPL_BLAS_LIBRARIES:PATH=/usr/local/libblas.a\n" \
      +"-DTPL_LAPACK_LIBRARIES:PATH=/usr/local/liblapack.a\n" \
      )

    writeStrToFile(testBaseDir+"/MPI_DEBUG_SS.config",
      "-DTPL_ENABLE_MPI:BOOL=ON\n" \
      +"-DTeuchos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON\n" \
      )

    modifiedFilesStr = ""

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      "--make-options=-j3 --ctest-options=-j5" \
      +" --default-builds=MPI_DEBUG --do-all --push " \
      +" --ss-extra-builds=MPI_DEBUG_SS --enable-packages=Teuchos,Phalanx" \
      ,
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsFinalPushPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "Phalanx of type SS is being excluded because it is not in the valid list of package types .PS.\n" \
      +"passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
      +"passed: Trilinos/MPI_DEBUG_SS: passed=100,notpassed=0\n" \
      +"0) MPI_DEBUG => passed: passed=100,notpassed=0\n" \
      +"2) MPI_DEBUG_SS => passed: passed=100,notpassed=0\n" \
      +"^DID PUSH\n" \
      )


  def test_ss_extra_builds_ex_only_fail(self):
    
    testName = "ss_extra_builds_ex_only_fail"

    testBaseDir = create_checkin_test_case_dir(testName, g_verbose)

    writeStrToFile(testBaseDir+"/COMMON.config",
      "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON\n" \
      +"-DBUILD_SHARED:BOOL=ON\n" \
      +"-DTPL_BLAS_LIBRARIES:PATH=/usr/local/libblas.a\n" \
      +"-DTPL_LAPACK_LIBRARIES:PATH=/usr/local/liblapack.a\n" \
      )

    writeStrToFile(testBaseDir+"/MPI_DEBUG_SS.config",
      "-DTPL_ENABLE_MPI:BOOL=ON\n" \
      +"-DTeuchos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON\n" \
      )

    modifiedFilesStr = ""

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      " --default-builds=MPI_DEBUG --send-email-to=" \
      +" --make-options=-j3 --ctest-options=-j5" \
      +" --do-all --push" \
      +" --ss-extra-builds=MPI_DEBUG_SS --enable-packages=ThyraCrazyStuff" \
      ,
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      ,
      \
      False,
      \
      "ThyraCrazyStuff of type EX is being excluded because it is not in the valid list of package types .PS.\n" \
      "ThyraCrazyStuff of type EX is being excluded because it is not in the valid list of package types .PS,SS.\n" \
      +"passed: Trilinos/MPI_DEBUG: skipped configure, build, test due to no enabled packages\n" \
      +"passed: Trilinos/MPI_DEBUG_SS: skipped configure, build, test due to no enabled packages\n" \
      +"There were no successfuly attempts to configure/build/test!\n" \
      +"  => A PUSH IS .NOT. READY TO BE PERFORMED!\n" \
      +"^PUSH FAILED\n" \
      +"^REQUESTED ACTIONS: FAILED\n" \
      )


  def test_ss_extra_builds_extra_builds_ex_only_pass(self):
    
    testName = "ss_extra_builds_extra_builds_ex_only_pass"

    testBaseDir = create_checkin_test_case_dir(testName, g_verbose)

    writeStrToFile(testBaseDir+"/COMMON.config",
      "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON\n" \
      +"-DBUILD_SHARED:BOOL=ON\n" \
      +"-DTPL_BLAS_LIBRARIES:PATH=/usr/local/libblas.a\n" \
      +"-DTPL_LAPACK_LIBRARIES:PATH=/usr/local/liblapack.a\n" \
      )

    writeStrToFile(testBaseDir+"/MPI_DEBUG_SS.config",
      "-DTPL_ENABLE_MPI:BOOL=ON\n" \
      +"-DTeuchos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON\n" \
      )

    writeStrToFile(testBaseDir+"/SERIAL_RELEASE_SS.config",
      "-DTPL_ENABLE_MPI:BOOL=OFF\n" \
      +"-DTeuchos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON\n" \
      )

    modifiedFilesStr = ""

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      " --default-builds=MPI_DEBUG" \
      +" --make-options=-j3 --ctest-options=-j5" \
      +" --do-all --push" \
      +" --ss-extra-builds=MPI_DEBUG_SS" \
      +" --extra-builds=SERIAL_RELEASE_SS" \
      +" --enable-packages=ThyraCrazyStuff" \
      ,
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsFinalPushPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "ThyraCrazyStuff of type EX is being excluded because it is not in the valid list of package types .PS.\n" \
      "ThyraCrazyStuff of type EX is being excluded because it is not in the valid list of package types .PS,SS.\n" \
      +"passed: Trilinos/MPI_DEBUG: skipped configure, build, test due to no enabled packages\n" \
      +"passed: Trilinos/MPI_DEBUG_SS: skipped configure, build, test due to no enabled packages\n" \
      +"passed: Trilinos/SERIAL_RELEASE_SS: passed=100,notpassed=0\n" \
      +"0) MPI_DEBUG => Skipped configure, build, test due to no enabled packages! => Does not affect push readiness!\n" \
      +"2) MPI_DEBUG_SS => Skipped configure, build, test due to no enabled packages! => Does not affect push readiness!\n" \
      +"3) SERIAL_RELEASE_SS => passed: passed=100,notpassed=0\n" \
      +"^DID PUSH\n" \
      )

  # E) Test intermediate states with rerunning to fill out


  # ToDo: Add test for pull followed by configure
  # ToDo: Add test for configure followed by build
  # ToDo: Add test for build followed by test


  def test_do_all_default_builds_mpi_debug_then_push_pass(self):

    testName = "do_all_default_builds_mpi_debug_then_push_pass"

    # Do the build/test only first (ready to push)
    g_test_do_all_default_builds_mpi_debug_pass(self, testName)

    # Do the push after the fact
    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      "--make-options=-j3 --ctest-options=-j5 --default-builds=MPI_DEBUG --push",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsFinalPushPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "0) MPI_DEBUG => passed: passed=100,notpassed=0\n" \
      "0) MPI_DEBUG Results:\n" \
      +"=> A PUSH IS READY TO BE PERFORMED!\n" \
      +"^DID PUSH: Trilinos:\n" \
      
      )


  def test_do_all_default_builds_mpi_debug_then_empty(self):

    testName = "do_all_default_builds_mpi_debug_then_push_pass"

    # Do the build/test only first (ready to push)
    g_test_do_all_default_builds_mpi_debug_pass(self, testName)

    # Check the status after (no action arguments)
    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      "--make-options=-j3 --ctest-options=-j5 --default-builds=MPI_DEBUG",
      \
      g_cmndinterceptsCurrentBranch \
      +"IT: eg diff --name-status origin/currentbranch; 0; 'eg diff passed'\n" 
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "0) MPI_DEBUG => passed: passed=100,notpassed=0\n" \
      "=> A PUSH IS READY TO BE PERFORMED!\n" \
      "^READY TO PUSH: Trilinos:\n"
      )


  def test_ss_extra_builds_ss_do_all_then_empty(self):

    testName = "test_ss_extra_builds_ss_do_all_then_empty"

    # --do-all without the push (ready to push)
    g_test_ss_extra_builds_ss_do_all_pass(self, testName)

    # Follow-up status (should be ready to push)

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      " --default-builds=MPI_DEBUG --ss-extra-builds=MPI_DEBUG_SS --enable-packages=Phalanx" \
      ,
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "0) MPI_DEBUG => passed: skipped configure, build, test due to no enabled packages\n" \
      +"2) MPI_DEBUG_SS => passed: passed=100,notpassed=0\n" \
      +"^READY TO PUSH\n" \
      )


  def test_ss_extra_builds_ss_do_all_then_push(self):

    testName = "test_ss_extra_builds_ss_do_all_then_push"

    # --do-all without the push (ready to push)
    g_test_ss_extra_builds_ss_do_all_pass(self, testName)

    # Follow-up status (should be ready to push)

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      " --default-builds=MPI_DEBUG --ss-extra-builds=MPI_DEBUG_SS --enable-packages=Phalanx" \
      +" --push"
      ,
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsFinalPushPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "0) MPI_DEBUG => passed: skipped configure, build, test due to no enabled packages\n" \
      +"2) MPI_DEBUG_SS => passed: passed=100,notpassed=0\n" \
      +"^DID PUSH\n" \
      )


  # ToDo: On all of these below check that the right files are being deleted!

  # ToDo: Add test for removing files on pull (fail immediately)
  # ToDo: Add test for removing files on configure (fail immediately)
  # ToDo: Add test for removing files on build (fail immediately)
  # ToDo: Add test for removing files on test (fail immediately)


  # F) Test various failing use cases


  def test_do_all_wrong_eg_version(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_wrong_eg_version",
      \
      "--do-all --eg-git-version-check" \
      ,
      \
      "IT: eg --version; 1; 'eg version wrong-version'\n" \
      ,
      \
      False,
      \
      "Error, the installed eg version wrong-version does not equal the official eg version "+g_officialEgVersion+"!\n" \
      ,
      egVersion=False
      )


  def test_wrong_eg_version_ignore(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_wrong_eg_vesion",
      \
      "--no-eg-git-version-check --skip-push-readiness-check" \
      ,
      \
      "IT: eg --version; 1; 'eg version wrong-version'\n" \
      +g_cmndinterceptsCurrentBranch \
      ,
      \
      True,
      \
      "WARNING: No actions were performed!\n" \
      "REQUESTED ACTIONS: PASSED\n" \
      ,
      egVersion=False
      )

  #NOTE: I would also like to check the git verion but I can't becuase my
  #command intercept system can't hanlde more than one line of output.


  def test_enable_packages_error(self):
    checkin_test_run_case(
      \
      self,
      \
      "enable_packages_error",
      \
      "--enable-packages=TEuchos" \
      ,
      \
      "" \
      ,
      \
      False,
      \
      "Error, invalid package name TEuchos in --enable-packages=TEuchos." \
      "  The valid package names include: .*Teuchos, .*\n" \
      )


  def test_disable_packages_error(self):
    checkin_test_run_case(
      \
      self,
      \
      "disable_packages_error",
      \
      "--disable-packages=TEuchos" \
      ,
      \
      "" \
      ,
      \
      False,
      \
      "Error, invalid package name TEuchos in --disable-packages=TEuchos." \
      "  The valid package names include: .*Teuchos, .*\n" \
      )


  def test_do_all_local_do_all(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_local_do_all",
      \
      "--do-all --local-do-all" \
      ,
      \
      "" \
      ,
      \
      False,
      \
      "Error, you can not use --do-all and --local-do-all together!\n" \
      )


  def test_do_all_allow_no_pull(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_allow_no_pull",
      \
      "--do-all --allow-no-pull" \
      ,
      \
      "" \
      ,
      \
      False,
      \
      "Error, you can not use --do-all and --allow-no-pull together!\n" \
      )


  def test_do_all_default_builds_mpi_debug_unstaged_changed_files_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_default_builds_mpi_debug_pull_fail",
      \
      "--do-all",
      g_cmndinterceptsCurrentBranch \
      +"IT: eg status; 0; 'Changed but not updated'\n" \
      ,
      \
      False,
      \
      "ERROR: There are changed unstaged uncommitted files => cannot continue!\n" \
      "Update failed!\n" \
      "Not running any build/test cases because the update (pull) failed!\n" \
      "A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      "INITIAL PULL FAILED: Trilinos:\n"
      )


  def test_do_all_default_builds_mpi_debug_staged_uncommitted_files_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_default_builds_mpi_debug_pull_fail",
      \
      "--do-all",
      g_cmndinterceptsCurrentBranch \
      +"IT: eg status; 0; 'Changes ready to be committed'\n" \
      ,
      \
      False,
      \
      "ERROR: There are changed staged uncommitted files => cannot continue!\n" \
      "Update failed!\n" \
      "Not running any build/test cases because the update (pull) failed!\n" \
      "A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      "INITIAL PULL FAILED: Trilinos:\n"
      )


  def test_do_all_default_builds_mpi_debug_unknown_files_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_default_builds_mpi_debug_pull_fail",
      \
      "--do-all",
      g_cmndinterceptsCurrentBranch \
      +"IT: eg status; 0; 'Newly created unknown files'\n" \
      ,
      \
      False,
      \
      "ERROR: There are newly created uncommitted files => Cannot continue!\n" \
      "Update failed!\n" \
      "Not running any build/test cases because the update (pull) failed!\n" \
      "A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      "INITIAL PULL FAILED: Trilinos:\n"
      )


  def test_do_all_default_builds_mpi_debug_pull_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_default_builds_mpi_debug_pull_fail",
      \
      "--do-all",
      g_cmndinterceptsCurrentBranch \
      +"IT: eg status; 0; '(on master branch)'\n" \
      +"IT: eg pull; 1; 'eg pull failed'\n" \
      ,
      \
      False,
      \
      "Pull failed!\n" \
      "Update failed!\n" \
      "Skipping getting list of modified files because pull failed!\n" \
      "A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      "INITIAL PULL FAILED: Trilinos:\n"
      )


  def test_illegal_enables_fail(self):
    
    testName = "illegal_enables_fail"

    testBaseDir = create_checkin_test_case_dir(testName, g_verbose)

    writeStrToFile(testBaseDir+"/COMMON.config",
      "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON\n" \
      +"-DBUILD_SHARED:BOOL=ON\n" \
      +"-DTPL_BLAS_LIBRARIES:PATH=/usr/local/libblas.a\n" \
      +"-DTPL_LAPACK_LIBRARIES:PATH=/usr/local/liblapack.a\n" \
      +"-DTPL_ENABLE_BOOST:BOOL=ON\n" \
      +"-DTrilinos_ENABLE_TriKota:BOOL=ON\n" \
      +"-DTrilinos_ENABLE_WebTrilinos=ON\n" \
      )

    writeStrToFile(testBaseDir+"/MPI_DEBUG.config",
      "-DTPL_ENABLE_MPI:BOOL=ON\n" \
      +"-DTPL_ENABLE_CUDA:BOOL=ON\n" \
      +"-DTrilinos_ENABLE_STK:BOOL=ON\n" \
      +"-DTrilinos_ENABLE_Phalanx=ON\n" \
      +"-DTrilinos_ENABLE_Sundance=OFF\n" \
      )

    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      "--default-builds=MPI_DEBUG --configure --allow-no-pull",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      "ERROR: Illegal TPL enable -DTPL_ENABLE_BOOST:BOOL=ON in ../COMMON.config!\n" \
      +"ERROR: Illegal enable -DTrilinos_ENABLE_TriKota:BOOL=ON in ../COMMON.config!\n" \
      +"ERROR: Illegal enable -DTrilinos_ENABLE_WebTrilinos=ON in ../COMMON.config!\n" \
      +"ERROR: Illegal TPL enable -DTPL_ENABLE_CUDA:BOOL=ON in ../MPI_DEBUG.config!\n" \
      +"ERROR: Illegal TPL enable -DTPL_ENABLE_MPI:BOOL=ON in ../MPI_DEBUG.config!\n" \
      +"ERROR: Illegal enable -DTrilinos_ENABLE_STK:BOOL=ON in ../MPI_DEBUG.config!\n" \
      +"ERROR: Illegal enable -DTrilinos_ENABLE_Phalanx=ON in ../MPI_DEBUG.config!\n" \
      +"SKIPPED: MPI_DEBUG configure skipped because pre-configure failed (see above)!\n" \
      +"0) MPI_DEBUG => FAILED: pre-configure failed => Not ready to push!\n" \
      +"Configure: FAILED\n" \
      +"FAILED CONFIGURE/BUILD/TEST: Trilinos:\n" \
      +"REQUESTED ACTIONS: FAILED\n" \
      ,
      failRegexStrList = \
      "ERROR: Illegal enable -DTrilinos_ENABLE_Sundance=OFF\n" # Package disables are okay but not great
      )

    self.assertEqual(os.path.exists(testBaseDir+"/MPI_DEBUG/do-configure.base"), False)
    self.assertEqual(os.path.exists(testBaseDir+"/MPI_DEBUG/do-configure"), False)


  def test_do_all_default_builds_mpi_debug_configure_fail(self):
    checkin_test_run_case(
      self,
      \
      "do_all_default_builds_mpi_debug_configure_fail",
      \
      "--do-all --default-builds=MPI_DEBUG",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +"IT: \./do-configure; 1; 'do-configure failed'\n" \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      g_expectedRegexUpdateWithBuildCasePasses \
      +"Configure failed returning 1!\n" \
      +"The configure FAILED!\n" \
      +"The build was never attempted!\n" \
      +"The tests where never even run!\n" \
      +"FAILED: configure failed\n" \
      +"A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      +"FAILED CONFIGURE/BUILD/TEST: Trilinos:\n" \
      )


  def test_do_all_default_builds_mpi_debug_build_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_default_builds_mpi_debug_build_fail",
      \
      "--do-all --default-builds=MPI_DEBUG --make-options=-j3 --ctest-options=-j5",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +"IT: \./do-configure; 0; 'do-configure passed'\n" \
      +"IT: make -j3; 1; 'make filed'\n" \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      g_expectedRegexUpdateWithBuildCasePasses \
      +g_expectedRegexConfigPasses \
      +g_expectedRegexTestNotRun \
      +g_expectedRegexBuildFailed \
      +"0) MPI_DEBUG => FAILED: build failed => Not ready to push!\n" \
      +"1) SERIAL_RELEASE => Test case SERIAL_RELEASE was not run! => Does not affect push readiness!\n" \
      +g_expectedCommonOptionsSummary \
      +"A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      +"^FAILED CONFIGURE/BUILD/TEST: Trilinos:\n"
      )


  def test_do_all_default_builds_mpi_debug_test_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_default_builds_mpi_debug_test_fail",
      \
      "--do-all --default-builds=MPI_DEBUG --make-options=-j3 --ctest-options=-j5",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +"IT: \./do-configure; 0; 'do-configure passed'\n" \
      +"IT: make -j3; 0; 'make passed'\n" \
      +"IT: ctest -j5; 1; '80% tests passed, 20 tests failed out of 100.\n" \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsSendFinalEmail \
      ,      \
      False,
      \
      g_expectedRegexUpdateWithBuildCasePasses \
      +g_expectedRegexConfigPasses \
      +g_expectedRegexBuildPasses \
      +"FAILED: ctest failed returning 1!\n" \
      +"testResultsLine = .80% tests passed, 20 tests failed out of 100.\n" \
      +"0) MPI_DEBUG => FAILED: passed=80,notpassed=20\n" \
      +"1) SERIAL_RELEASE => Test case SERIAL_RELEASE was not run! => Does not affect push readiness!\n" \
      +g_expectedCommonOptionsSummary \
      +"Test: FAILED\n" \
      +"A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      +"^FAILED CONFIGURE/BUILD/TEST: Trilinos:\n" \
      )


  def test_do_all_push_default_builds_mpi_debug_final_pull_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_push_default_builds_mpi_debug_final_pull_fail",
      \
      "--default-builds=MPI_DEBUG --make-options=-j3 --ctest-options=-j5" \
      " --do-all --push" \
      ,
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +"IT: eg pull && eg rebase --against origin/currentbranch; 1; 'final eg pull FAILED'\n" \
      +"IT: eg log --oneline currentbranch \^origin/currentbranch; 0; '54321 Only one commit'\n" \
      +g_cmndinterceptsSendFinalEmail \
      ,      \
      False,
      \
      g_expectedRegexUpdateWithBuildCasePasses \
      +g_expectedRegexConfigPasses \
      +g_expectedRegexBuildPasses \
      +g_expectedRegexTestPasses \
      +g_expectedCommonOptionsSummary \
      +"A PUSH IS READY TO BE PERFORMED!\n" \
      +"Final update failed!\n" \
      +"Skippng appending test results due to prior errors!\n" \
      +"Not performing push due to prior errors!\n" \
      +"FINAL PULL FAILED: Trilinos:\n" \
      )


  def test_do_all_push_default_builds_mpi_debug_final_commit_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_push_default_builds_mpi_debug_final_commit_fail",
      \
      "--default-builds=MPI_DEBUG --make-options=-j3 --ctest-options=-j5" \
      " --do-all --push" \
      ,
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +"IT: eg pull && eg rebase --against origin/currentbranch; 0; 'final eg pull and rebase passed'\n" \
      +g_cmnginterceptsEgLogCmnds \
      +"IT: eg commit --amend -F .*; 1; 'Amending the last commit FAILED'\n" \
      +"IT: eg log --oneline currentbranch \^origin/currentbranch; 0; '54321 Only one commit'\n" \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      g_expectedRegexUpdateWithBuildCasePasses \
      +g_expectedRegexConfigPasses \
      +g_expectedRegexBuildPasses \
      +g_expectedRegexTestPasses \
      +g_expectedCommonOptionsSummary \
      +"A PUSH IS READY TO BE PERFORMED!\n" \
      +"Final update passed!\n" \
      +"Attempting to amend the final commmit message ...\n" \
      +"Appending test results to last commit failed!\n" \
      +"Not performing push due to prior errors!\n" \
      +"AMEND COMMIT FAILED: Trilinos:\n" \
      )


  def test_do_all_push_default_builds_mpi_debug_push_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_push_default_builds_mpi_debug_push_fail",
      \
      "--default-builds=MPI_DEBUG --make-options=-j3 --ctest-options=-j5" \
      " --do-all --push" \
      ,
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +"IT: eg pull && eg rebase --against origin/currentbranch; 0; 'final eg pull and rebase passed'\n" \
      +g_cmnginterceptsEgLogCmnds \
      +"IT: eg commit --amend -F .*; 0; 'Amending the last commit passed'\n" \
      +"IT: eg log --oneline currentbranch \^origin/currentbranch; 0; '54321 Only one commit'\n" \
      +"IT: cat modifiedFiles.out; 0; 'M\tpackages/teuchos/CMakeLists.txt'\n"\
      +"IT: eg push; 1; 'push FAILED'\n"
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      g_expectedRegexUpdateWithBuildCasePasses \
      +g_expectedRegexConfigPasses \
      +g_expectedRegexBuildPasses \
      +g_expectedRegexTestPasses \
      +g_expectedCommonOptionsSummary \
      +"A PUSH IS READY TO BE PERFORMED!\n" \
      +"Final update passed!\n" \
      +"Appending test results to last commit passed!\n" \
      +"Push failed!\n" \
      +"PUSH FAILED: Trilinos:\n" \
      )


  def test_do_all_push_no_local_commits_push_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_push_no_local_commits_push_fail",
      \
      "--make-options=-j3 --ctest-options=-j5 --do-all --push",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +"IT: eg pull && eg rebase --against origin/currentbranch; 0; 'final eg pull and rebase passed'\n" \
      +"IT: eg cat-file -p HEAD; 0; 'This is the last commit message'\n" \
      +"IT: eg log --oneline currentbranch \^origin/currentbranch; 0; ''\n" \
      +"IT: eg log --pretty=format:'%h' currentbranch\^ \^origin/currentbranch; 0; ''\n" \
      +"IT: eg log --oneline currentbranch \^origin/currentbranch; 0; '54321 Only one commit'\n" \
      +"IT: cat modifiedFiles.out; 0; ''\n"\
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      g_expectedRegexUpdateWithBuildCasePasses \
      +g_expectedRegexConfigPasses \
      +g_expectedRegexBuildPasses \
      +g_expectedRegexTestPasses \
      +"0) MPI_DEBUG => passed: passed=100,notpassed=0\n" \
      +"1) SERIAL_RELEASE => passed: passed=100,notpassed=0\n" \
      +"=> A PUSH IS READY TO BE PERFORMED!\n" \
      +"No local commits exit!\n" \
      +"Skipping amending last commit because there are no local commits!\n" \
      +"Attempting to do the push ...\n" \
      +"Skipping push to .. because there are no changes!\n" \
      +"Push failed because the push was never attempted!\n" \
      +"^PUSH FAILED: Trilinos:\n" \
      )


  def test_do_all_default_builds_mpi_debug_push_no_tests_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_default_builds_mpi_debug_push_no_tests_fail",
      \
      "--make-options=-j3 --ctest-options=-j5 --default-builds=MPI_DEBUG --do-all --push",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildPasses \
      +"IT: ctest -j5; 0; 'No tests were found!!!'\n" \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      g_expectedRegexConfigPasses \
      +g_expectedRegexBuildPasses \
      +"No tests failed!\n"\
      +"CTest was invoked but no tests were run!\n"\
      +"At least one of the actions (update, configure, built, test) failed or was not performed correctly!\n" \
       +"0) MPI_DEBUG => FAILED: no tests run\n" \
      +"=> A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      +"^FAILED CONFIGURE/BUILD/TEST: Trilinos:\n" \
      +"REQUESTED ACTIONS: FAILED\n" \
      )


  def test_local_do_all_default_builds_mpi_debug_push_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "local_do_all_default_builds_mpi_debug_push_fail",
      \
      "--make-options=-j3 --ctest-options=-j5 --default-builds=MPI_DEBUG --local-do-all --push",
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      g_expectedRegexConfigPasses \
      +g_expectedRegexBuildPasses \
      +g_expectedRegexTestPasses \
      +"0) MPI_DEBUG => passed: passed=100,notpassed=0\n" \
      +"A current successful pull does \*not\* exist => Not ready for final push!\n"\
      +"=> A PUSH IS \*NOT\* READY TO BE PERFORMED!\n"\
      +"^ABORTED COMMIT/PUSH: Trilinos:\n" \
      +"REQUESTED ACTIONS: FAILED\n" \
      )


  def test_extra_repo_1_no_changes_do_all_push_fail(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_no_changes_do_all_push_fail",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --default-builds=MPI_DEBUG --do-all --push", \
      \
      g_cmndinterceptsExtraRepo1DoAllUpToPush \
      +g_cmndinterceptsCatModifiedFilesNoChanges \
      +g_cmndinterceptsCatModifiedFilesPreCoprightTrilinosNoChanges \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      "Skipping push to .. because there are no changes!\n" \
      "Skipping push to .preCopyrightTrilinos. because there are no changes!\n" \
      +"Push failed because the push was never attempted!\n" \
      +"PUSH FAILED: Trilinos:\n" \
      +"REQUESTED ACTIONS: FAILED\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_send_email_only_on_failure_do_all_mpi_debug_build_configure_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "send_email_only_on_failure_do_all_mpi_debug_build_configure_fail",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      +" --send-email-only-on-failure" \
      +" --do-all --push" \
      ,
      \
      g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsPullPasses \
      +"IT: \./do-configure; 1; 'do-configure failed'\n" \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      g_expectedRegexUpdateWithBuildCasePasses \
      +"Configure failed returning 1!\n" \
      +"The configure FAILED!\n" \
      +"The build was never attempted!\n" \
      +"The tests where never even run!\n" \
      +"FAILED: configure failed\n" \
      +g_expectedRegexConfigPasses \
      +g_expectedRegexBuildPasses \
      +g_expectedRegexTestPasses \
      +g_expectedCommonOptionsSummary \
      +"Running: mailx -s .FAILED: Trilinos/MPI_DEBUG: configure failed. bogous@somwhere.com\n" \
      +"SERIAL_RELEASE: Skipping sending build/test case email because it passed and --send-email-only-on-failure was set!\n" \
      +"Running: mailx -s .FAILED CONFIGURE/BUILD/TEST: Trilinos: .* bogous@somwhere.com\n" \
      )


  def test_extra_repo_1_mispell_repo_fail(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "test_extra_repo_1_mispell_repo_fail",
      \
      " --extra-repos=preCopyrightTrilinosMispell", \
      \
      g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      "Error, the specified git repo .preCopyrightTrilinosMispell. directory .*preCopyrightTrilinosMispell. does not exist!\n"
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_initial_trilinos_pull_fail(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_initial_trilinos_pull_fail",
      \
      " --extra-repos=preCopyrightTrilinos --default-builds=MPI_DEBUG --pull", \
      \
      "IT: .*cmake .+ -P .+/TribitsDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
      +g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsPullOnlyFails \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      "pullInitial.out\n" \
      +"Pull failed!\n" \
      +"Skipping getting list of modified files because pull failed!\n" \
      +"INITIAL PULL FAILED: Trilinos:\n" \
      +"REQUESTED ACTIONS: FAILED\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_initial_extra_repo_pull_fail(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_initial_extra_repo_pull_fail",
      \
      " --extra-repos=preCopyrightTrilinos --default-builds=MPI_DEBUG --pull", \
      \
      "IT: .*cmake .+ -P .+/TribitsDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
      +g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyFails \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      "pullInitial.out\n" \
      "pullInitial.preCopyrightTrilinos.out\n" \
      +"Pull failed!\n" \
      +"Skipping getting list of modified files because pull failed!\n" \
      +"INITIAL PULL FAILED: Trilinos:\n" \
      +"REQUESTED ACTIONS: FAILED\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_extra_pull_trilinos_fail(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_extra_pull_trilinos_fail",
      \
      " --extra-repos=preCopyrightTrilinos --default-builds=MPI_DEBUG --pull --extra-pull-from=ssg:master", \
      \
      "IT: .*cmake .+ -P .+/TribitsDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
      +g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyFails \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      "pullInitial.out\n" \
      "pullInitial.preCopyrightTrilinos.out\n" \
      "pullInitialExtra.out\n" \
      +"Pull failed!\n" \
      +"Skipping getting list of modified files because pull failed!\n" \
      +"INITIAL PULL FAILED: Trilinos:\n" \
      +"REQUESTED ACTIONS: FAILED\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_extra_pull_extra_repo_fail(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_extra_pull_extra_repo_fail",
      \
      " --extra-repos=preCopyrightTrilinos --default-builds=MPI_DEBUG --pull --extra-pull-from=ssg:master", \
      \
      "IT: .*cmake .+ -P .+/TribitsDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
      +g_cmndinterceptsCurrentBranch \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsStatusPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyPasses \
      +g_cmndinterceptsPullOnlyFails \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      "pullInitial.out\n" \
      "pullInitial.preCopyrightTrilinos.out\n" \
      "pullInitialExtra.out\n" \
      "pullInitialExtra.preCopyrightTrilinos.out\n" \
      "Pull failed!\n" \
      "Skipping getting list of modified files because pull failed!\n" \
      "INITIAL PULL FAILED: Trilinos:\n" \
      "REQUESTED ACTIONS: FAILED\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_do_all_final_pull_trilinos_fails(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_do_all_final_pull_trilinos_fails",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --default-builds=MPI_DEBUG --do-all --push", \
      \
      g_cmndinterceptsExtraRepo1DoAllThroughTest \
      +g_cmndinterceptsFinalPullRebaseFails \
      +g_cmndinterceptsLogCommitsPasses \
      +g_cmndinterceptsLogCommitsPasses \
      +g_cmndinterceptsSendFinalEmail
      ,
      \
      False,
      \
      "pullFinal.out\n" \
      "FINAL PULL FAILED: Trilinos:\n" \
      "REQUESTED ACTIONS: FAILED\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_do_all_final_pull_extra_repo_fails(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_do_all_final_pull_extra_repo_fails",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --default-builds=MPI_DEBUG --do-all --push", \
      \
      g_cmndinterceptsExtraRepo1DoAllThroughTest \
      +g_cmndinterceptsFinalPullRebasePasses \
      +g_cmndinterceptsFinalPullRebaseFails \
      +g_cmndinterceptsLogCommitsPasses \
      +g_cmndinterceptsLogCommitsPasses \
      +g_cmndinterceptsSendFinalEmail
      ,
      \
      False,
      \
      "pullFinal.out\n" \
      "pullFinal.preCopyrightTrilinos.out\n" \
      "Final update failed!\n" \
      "FINAL PULL FAILED: Trilinos:\n" \
      "REQUESTED ACTIONS: FAILED\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_do_all_final_amend_trilinos_fails(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_do_all_final_amend_trilinos_fails",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --default-builds=MPI_DEBUG --do-all --push", \
      \
      g_cmndinterceptsExtraRepo1DoAllThroughTest \
      +g_cmndinterceptsFinalPullRebasePasses \
      +g_cmndinterceptsFinalPullRebasePasses \
      +g_cmnginterceptsEgLogCmnds \
      +g_cmndinterceptsAmendCommitFails \
      +g_cmndinterceptsLogCommitsPasses \
      +g_cmndinterceptsLogCommitsPasses \
      +g_cmndinterceptsSendFinalEmail
      ,
      \
      False,
      \
      "commitFinalBody.out\n" \
      "Appending test results to last commit failed!\n" \
      "AMEND COMMIT FAILED: Trilinos:\n" \
      "REQUESTED ACTIONS: FAILED\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_do_all_final_amend_extra_repo_fails(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_do_all_final_amend_extra_repo_fails",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --default-builds=MPI_DEBUG --do-all --push", \
      \
      g_cmndinterceptsExtraRepo1DoAllThroughTest \
      +g_cmndinterceptsFinalPullRebasePasses \
      +g_cmndinterceptsFinalPullRebasePasses \
      +g_cmnginterceptsEgLogCmnds \
      +g_cmndinterceptsAmendCommitPasses \
      +g_cmnginterceptsEgLogCmnds \
      +g_cmndinterceptsAmendCommitFails \
      +g_cmndinterceptsLogCommitsPasses \
      +g_cmndinterceptsLogCommitsPasses \
      +g_cmndinterceptsSendFinalEmail
      ,
      \
      False,
      \
      "commitFinalBody.out\n" \
      "commitFinalBody.preCopyrightTrilinos.out\n" \
      "Appending test results to last commit failed!\n" \
      "AMEND COMMIT FAILED: Trilinos:\n" \
      "REQUESTED ACTIONS: FAILED\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_do_all_final_push_trilinos_fails(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_do_all_final_push_trilinos_fails",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --default-builds=MPI_DEBUG --do-all --push", \
      \
      g_cmndinterceptsExtraRepo1DoAllUpToPush \
      +g_cmndinterceptsCatModifiedFilesPasses \
      +g_cmndinterceptsPushOnlyFails \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      "Push failed!\n" \
      "PUSH FAILED: Trilinos:\n" \
      "REQUESTED ACTIONS: FAILED\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


  def test_extra_repo_1_do_all_final_push_extra_repo_fails(self):
    projectDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_do_all_final_push_trilinos_fails",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --default-builds=MPI_DEBUG --do-all --push", \
      \
      g_cmndinterceptsExtraRepo1DoAllUpToPush \
      +g_cmndinterceptsCatModifiedFilesPasses \
      +g_cmndinterceptsPushOnlyPasses \
      +g_cmndinterceptsCatModifiedFilesPreCoprightTrilinosPasses \
      +g_cmndinterceptsPushOnlyFails \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      "Push failed!\n" \
      "PUSH FAILED: Trilinos:\n" \
      "REQUESTED ACTIONS: FAILED\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+projectDepsXmlFileOverride ]
      )


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(testCheckinTest))
    return suite


from optparse import OptionParser


if __name__ == '__main__':

  # Look for --with-cmake=??? argument and process and remove it
  global g_withCmake
  g_withCmake = "cmake"
  args = []
  for arg in sys.argv:
    arg_find_cmake = arg.find("--with-cmake")
    if arg_find_cmake == 0:
      g_withCmake = arg.split("=")[1]
    else:
      args.append(arg)
  sys.argv = args
  
  if os.path.exists(g_checkin_test_tests_dir):
    echoRunSysCmnd("rm -rf "+g_checkin_test_tests_dir, verbose=g_verbose)

  unittest.main()
