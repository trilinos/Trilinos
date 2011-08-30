#######################################
# Unit testing code for CheckinTest.py #
########################################


from CheckinTest import *
import unittest


class MockOptions:
  def __init__(self):
    self.enableAllPackages = 'auto'


scriptsDir = getScriptBaseDir()


#
# Test formatMinutesStr
#


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


#
# Test formatMinutesStr
#


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


#
# Test extractPackageEnablesFromChangeStatus
#


trilinosDepsXmlFileDefaultOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.gold.xml"
trilinosDependenciesDefault = getTrilinosDependenciesFromXmlFile(trilinosDepsXmlFileDefaultOverride)


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

    extractPackageEnablesFromChangeStatus(updateOutputStr, options, "",
      enablePackagesList, False, trilinosDependenciesDefault)

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

    extractPackageEnablesFromChangeStatus(updateOutputStr, options, "",
      enablePackagesList, False, trilinosDependenciesDefault)

    self.assertEqual( options.enableAllPackages, 'auto' )
    self.assertEqual( enablePackagesList, [u'TrilinosFramework', u'Stratimikos', u'ThyraCoreLibs', u'Tpetra'] )


  def test_extra_repo(self):

    updateOutputStr = """
M	ExtraTrilinosPackages.cmake
M	stalix/README
"""
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    trilinosDependenciesLocal = getTrilinosDependenciesFromXmlFile(trilinosDepsXmlFileOverride)

    options = MockOptions()
    enablePackagesList = []

    extractPackageEnablesFromChangeStatus(updateOutputStr, options, "preCopyrightTrilinos",
      enablePackagesList, False, trilinosDependenciesLocal)

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


#############################################################################
#
#                     Test checkin-test.py script
#
#############################################################################


# Test Data

g_cmndinterceptsCurrentBranch = \
  "IT: eg branch; 0; '* currentbranch'\n"

g_cmndinterceptsStatusPasses = \
  "IT: eg status; 0; '(on master branch)'\n"

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
  "IT: cmake .+ -P .+/PackageArchDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
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


def assertGrepFileForRegexStrList(testObject, testName, fileName, regexStrList, verbose):
  assert(os.path.isfile(fileName))
  for regexToFind in regexStrList.strip().split('\n'):
    if regexToFind == "": continue
    foundRegex = getCmndOutput("grep '"+regexToFind+"' "+fileName, True, False)
    if verbose or not foundRegex:
      print "\ncheckin_test::"+testName+": In '"+fileName+"' look for regex '"+regexToFind+"' ...", 
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
      print "\ncheckin_test::"+testName+": In '"+fileName \
        +"' assert not exist regex '"+regexToFind+"' ... '"+foundRegex+"'", 
      if foundRegex: print ": FAILED"
      else: print ": PASSED"
    testObject.assertEqual(foundRegex, "")


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

    
    cmnd = scriptsDir + "/../../checkin-test.py" \
     +" --no-eg-git-version-check" \
     +" --trilinos-src-dir="+scriptsDir+"/../DependencyUnitTests/MockTrilinos" \
     +" --send-email-to=bogous@somwhere.com" \
     + " " + optionsStr
    # NOTE: Above, we want to turn off the eg/git version tests since we want
    # these unit tests to run on machines that do not have the official
    # versions (e.g. the SCICO LAN) but where the versions might be okay.
    # Also, we have to point to the static mock Trilinos source directory
    # also so that preCopyrighTrilinos will show up as an extra repo.
    
    # C) Set up the command intercept file

    baseCmndInterceptsStr = \
      "FT: .*checkin-test-impl\.py.*\n" \
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
def g_test_do_all_without_serial_release_pass(testObject, testName):
  checkin_test_run_case(
    \
    testObject,
    \
    testName,
    \
    "--make-options=-j3 --ctest-options=-j5 --without-serial-release --do-all",
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
  fileFailRegexStrList=[], modifiedFilesStr="" \
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
    ,
    filePassRegexStrList
    ,
    fileFailRegexStrList=fileFailRegexStrList
    )


def checkin_test_configure_enables_test(testObject, testName, optionsStr, regexListStr, \
  notRegexListStr="", modifiedFilesStr="" \
  ):
  checkin_test_configure_test(
     testObject,
     testName,
     "--without-serial-release "+optionsStr,
     [("MPI_DEBUG/do-configure", regexListStr)],
     [("MPI_DEBUG/do-configure", notRegexListStr)],
     modifiedFilesStr,
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
       "\-DTPL_ENABLE_Pthread:BOOL=OFF\n" \
       +"\-DTPL_ENABLE_BinUtils:BOOL=OFF\n" \
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
       "\-DTPL_ENABLE_Pthread:BOOL=OFF\n" \
       +"\-DTPL_ENABLE_BinUtils:BOOL=OFF\n" \
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
    eg = scriptsDir+"/../DependencyUnitTests/MockTrilinos/commonTools/git/eg"
    checkin_test_run_case(
      \
      self,
      \
      "do_all_no_eg_installed",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      +" --without-serial-release" \
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
      "Setting to default eg in source tree '.*/commonTools/git/eg'\n" \
      ,
      inPathEg=False, egVersion=False
      )


  def test_do_all_without_serial_release_pass(self):
    g_test_do_all_without_serial_release_pass(self, "do_all_without_serial_release_pass")


  def test_local_do_all_without_serial_release_pass(self):
    checkin_test_run_case(
      \
      self,
      \
      "local_do_all_without_serial_release_pass",
      \
      "--make-options=-j3 --ctest-options=-j5 --without-serial-release" \
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


  def test_do_all_without_serial_release_test_fail_force_push_pass(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_without_serial_release_test_fail_force_push_pass",
      \
      "--make-options=-j3 --ctest-options=-j5 --without-serial-release" \
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


  def test_do_all_without_serial_release_then_wipe_clean_pull_pass(self):

    testName = "do_all_without_serial_release_then_from_scratch_pull_pass"

    # Do the build/test only first (ready to push)
    g_test_do_all_without_serial_release_pass(self, testName)

    # Do the push after the fact
    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      "--make-options=-j3 --ctest-options=-j5 --without-serial-release" \
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
      "Skipping configure because no packages are enabled and --abort-gracefully-if-no-enables!\n" \
      +"subjectLine = .passed: Trilinos/MPI_DEBUG: skipped configure, build, test due to no enabled packages.\n" \
      +"subjectLine = .passed: Trilinos/SERIAL_RELEASE: skipped configure, build, test due to no enabled packages.\n" \
      +"0) MPI_DEBUG => passed: skipped configure, build, test due to no enabled packages => Not ready to push!\n" \
      +"1) SERIAL_RELEASE => passed: skipped configure, build, test due to no enabled packages => Not ready to push!\n" \
      +"MPI_DEBUG: Skipping sending build/test case email because there were no enables and --abort-gracefully-if-no-enables was set!\n"
      +"SERIAL_RELEASE: Skipping sending build/test case email because there were no enables and --abort-gracefully-if-no-enables was set!\n"
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
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_configure_pass",
      \
      "--extra-repos=preCopyrightTrilinos --allow-no-pull --without-serial-release --enable-packages=Stalix --configure", \
      \
      "IT: cmake .+ -P .+/PackageArchDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
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
      +"trilinosDepsXmlFileOverride="+trilinosDepsXmlFileOverride+"\n" \
      +"Enabling only the explicitly specified packages .Stalix. ...\n" \
      +"Trilinos_ENABLE_Stalix:BOOL=ON\n" \
      +"Trilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos\n" \
      +"Enabled Packages: Stalix\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_implicit_enable_configure_pass(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_configure_pass",
      \
      "--extra-repos=preCopyrightTrilinos --allow-no-pull --without-serial-release --configure", \
      \
      "IT: cmake .+ -P .+/PackageArchDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
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
      +"trilinosDepsXmlFileOverride="+trilinosDepsXmlFileOverride+"\n" \
      +"Modified file: .packages/../preCopyrightTrilinos/teko/CMakeLists.txt.\n" \
      +"  => Enabling .Teko.!\n" \
      +"Trilinos_ENABLE_Teko:BOOL=ON\n" \
      +"Trilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos\n" \
      +"Enabled Packages: Teuchos, Teko\n" \
      ,
      \
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_do_all_push_pass(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_do_all_push_pass",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --without-serial-release --do-all --push", \
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
      +"Trilinos_ENABLE_Teko:BOOL=ON\n" \
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_pull_extra_pull_pass(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_pull_extra_pull_pass",
      \
      "--extra-repos=preCopyrightTrilinos --pull --extra-pull-from=somemachine:someotherbranch", \
      \
      "IT: cmake .+ -P .+/PackageArchDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_trilinos_changes_do_all_push_pass(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_trilinos_changes_do_all_push_pass",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --without-serial-release --do-all --push", \
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_extra_repo_changes_do_all_push_pass(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_extra_repo_changes_do_all_push_pass",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --without-serial-release --do-all --push", \
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_abort_gracefully_if_no_updates_no_updates_passes(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_extra_pull_abort_gracefully_if_no_updates_no_updates_passes(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_extra_pull_abort_gracefully_if_no_updates_main_repo_update(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_extra_pull_abort_gracefully_if_no_updates_extra_repo_update(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_extra_pull_abort_gracefully_if_no_updates_main_repo_extra_update(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_extra_pull_abort_gracefully_if_no_updates_extra_repo_extra_update(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
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
      "--without-serial-release",
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
      "--disable-packages=Tpetra,Sundance",
      \
      "\-DTrilinos_ENABLE_Teuchos:BOOL=ON\n" \
      +"\-DTrilinos_ENABLE_Tpetra:BOOL=OFF\n" \
      +"\-DTrilinos_ENABLE_Sundance:BOOL=OFF\n" \
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
      +" --disable-packages=Tpetra,Sundance",
      \
      "\-DTrilinos_ENABLE_TrilinosFramework:BOOL=ON\n" \
      +"\-DTrilinos_ENABLE_RTOp:BOOL=ON\n" \
      +"\-DTrilinos_ENABLE_Thyra:BOOL=ON\n" \
      +"\-DTrilinos_ENABLE_Tpetra:BOOL=OFF\n" \
      +"\-DTrilinos_ENABLE_Sundance:BOOL=OFF\n" \
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
      "\-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON\n" \
      +"\-DTrilinos_ENABLE_TrilinosFramework:BOOL=ON\n",
      modifiedFilesStr="M\tcmake/utils/AppendSet.cmake",
      )


  def test_enable_all_packages_on(self):
    checkin_test_configure_enables_test(
      self,
      "enable_all_packages_on",
      "--enable-all-packages=on",
      "\-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON\n" \
      +"\-DTrilinos_ENABLE_Teuchos:BOOL=ON\n",
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


  def test_without_serial_release_pull_only(self):
    checkin_test_run_case(
      self,
      \
      "without_serial_release_pull_only",
      \
      "--without-serial-release --pull",
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


  def test_without_serial_release_pull_skip_push_readiness_check(self):
    checkin_test_run_case(
      self,
      \
      "without_serial_release_pull_skip_push_readiness_check",
      \
      "--without-serial-release --pull --skip-push-readiness-check",
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


  def test_without_serial_release_pull_extra_pull_only(self):
    checkin_test_run_case(
      self,
      \
      "without_serial_release_pull_extra_pull_only",
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


  def test_without_serial_release_extra_pull_only(self):
    checkin_test_run_case(
      self,
      \
      "without_serial_release_extra_pull_only",
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


  def test_without_serial_release_configure_only(self):
    checkin_test_run_case(
      self,
      \
      "without_serial_release_configure_only",
      \
      "--without-serial-release --pull --configure",
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


  def test_without_serial_release_build_only(self):
    checkin_test_run_case(
      self,
      \
      "without_serial_release_build_only",
      \
      "--make-options=-j3 --without-serial-release --pull --configure --build",
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


  # D) Test --extra-builds


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


  # E) Test intermediate states with rerunning to fill out


  # ToDo: Add test for pull followed by configure
  # ToDo: Add test for configure followed by build
  # ToDo: Add test for build followed by test


  def test_do_all_without_serial_release_then_push_pass(self):

    testName = "do_all_without_serial_release_then_push_pass"

    # Do the build/test only first (ready to push)
    g_test_do_all_without_serial_release_pass(self, testName)

    # Do the push after the fact
    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      "--make-options=-j3 --ctest-options=-j5 --without-serial-release --push",
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


  def test_do_all_without_serial_release_then_push_pass(self):

    testName = "do_all_without_serial_release_then_push_pass"

    # Do the build/test only first (ready to push)
    g_test_do_all_without_serial_release_pass(self, testName)

    # Do the push after the fact
    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      "--make-options=-j3 --ctest-options=-j5 --without-serial-release --push" \
      +" --extra-pull-from=dummy:master" \
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
      "0) MPI_DEBUG => passed: passed=100,notpassed=0\n" \
      +"=> A PUSH IS READY TO BE PERFORMED!\n" \
      +"^DID PUSH: Trilinos:\n" \
      ,
      failRegexStrList = \
      "eg pull dummy master\n" \
      )


  def test_do_all_without_serial_release_then_empty(self):

    testName = "do_all_without_serial_release_then_push_pass"

    # Do the build/test only first (ready to push)
    g_test_do_all_without_serial_release_pass(self, testName)

    # Check the status after (no action arguments)
    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      "--make-options=-j3 --ctest-options=-j5 --without-serial-release",
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


  def test_do_all_without_serial_release_unstaged_changed_files_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_without_serial_release_pull_fail",
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


  def test_do_all_without_serial_release_staged_uncommitted_files_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_without_serial_release_pull_fail",
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


  def test_do_all_without_serial_release_unknown_files_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_without_serial_release_pull_fail",
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


  def test_do_all_without_serial_release_pull_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_without_serial_release_pull_fail",
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
      "--without-serial-release --configure --allow-no-pull",
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
      +"Skipping configure because pre-configure failed (see above)!\n" \
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


  def test_do_all_without_serial_release_configure_fail(self):
    checkin_test_run_case(
      self,
      \
      "do_all_without_serial_release_configure_fail",
      \
      "--do-all --without-serial-release",
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


  def test_do_all_without_serial_release_build_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_without_serial_release_build_fail",
      \
      "--do-all --without-serial-release --make-options=-j3 --ctest-options=-j5",
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


  def test_do_all_without_serial_release_test_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_without_serial_release_test_fail",
      \
      "--do-all --without-serial-release --make-options=-j3 --ctest-options=-j5",
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


  def test_do_all_push_without_serial_release_final_pull_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_push_without_serial_release_final_pull_fail",
      \
      "--without-serial-release --make-options=-j3 --ctest-options=-j5" \
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


  def test_do_all_push_without_serial_release_final_commit_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_push_without_serial_release_final_commit_fail",
      \
      "--without-serial-release --make-options=-j3 --ctest-options=-j5" \
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


  def test_do_all_push_without_serial_release_push_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_push_without_serial_release_push_fail",
      \
      "--without-serial-release --make-options=-j3 --ctest-options=-j5" \
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


  def test_do_all_without_serial_release_push_no_tests_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_without_serial_release_push_no_tests_fail",
      \
      "--make-options=-j3 --ctest-options=-j5 --without-serial-release --do-all --push",
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


  def test_local_do_all_without_serial_release_push_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "local_do_all_without_serial_release_push_fail",
      \
      "--make-options=-j3 --ctest-options=-j5 --without-serial-release --local-do-all --push",
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
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_no_changes_do_all_push_fail",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --without-serial-release --do-all --push", \
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
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
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_initial_trilinos_pull_fail(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_initial_trilinos_pull_fail",
      \
      " --extra-repos=preCopyrightTrilinos --without-serial-release --pull", \
      \
      "IT: cmake .+ -P .+/PackageArchDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_initial_extra_repo_pull_fail(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_initial_extra_repo_pull_fail",
      \
      " --extra-repos=preCopyrightTrilinos --without-serial-release --pull", \
      \
      "IT: cmake .+ -P .+/PackageArchDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_extra_pull_trilinos_fail(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_extra_pull_trilinos_fail",
      \
      " --extra-repos=preCopyrightTrilinos --without-serial-release --pull --extra-pull-from=ssg:master", \
      \
      "IT: cmake .+ -P .+/PackageArchDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_extra_pull_extra_repo_fail(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_extra_pull_extra_repo_fail",
      \
      " --extra-repos=preCopyrightTrilinos --without-serial-release --pull --extra-pull-from=ssg:master", \
      \
      "IT: cmake .+ -P .+/PackageArchDumpDepsXmlScript.cmake; 0; 'dump XML file passed'\n" \
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_do_all_final_pull_trilinos_fails(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_do_all_final_pull_trilinos_fails",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --without-serial-release --do-all --push", \
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_do_all_final_pull_extra_repo_fails(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_do_all_final_pull_extra_repo_fails",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --without-serial-release --do-all --push", \
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_do_all_final_amend_trilinos_fails(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_do_all_final_amend_trilinos_fails",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --without-serial-release --do-all --push", \
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_do_all_final_amend_extra_repo_fails(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_do_all_final_amend_extra_repo_fails",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --without-serial-release --do-all --push", \
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_do_all_final_push_trilinos_fails(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_do_all_final_push_trilinos_fails",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --without-serial-release --do-all --push", \
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


  def test_extra_repo_1_do_all_final_push_extra_repo_fails(self):
    trilinosDepsXmlFileOverride=scriptsDir+"/UnitTests/TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml"
    checkin_test_run_case(
      \
      self,
      \
      "extra_repo_1_do_all_final_push_trilinos_fails",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --extra-repos=preCopyrightTrilinos --without-serial-release --do-all --push", \
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
      envVars = [ "CHECKIN_TEST_DEPS_XML_FILE_OVERRIDE="+trilinosDepsXmlFileOverride ]
      )


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(testCheckinTest))
    return suite


if __name__ == '__main__':
  if os.path.exists(g_checkin_test_tests_dir):
    echoRunSysCmnd("rm -rf "+g_checkin_test_tests_dir, verbose=g_verbose)
  unittest.main()
