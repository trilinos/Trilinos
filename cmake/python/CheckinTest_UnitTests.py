########################################
# Unit testing code for CheckinTest.py #
########################################


from CheckinTest import *
import unittest


class MockOptions:
  def __init__(self):
    self.enableAllPackages = 'default'


#
# Test isGlobalBuildFile
#

class test_isGlobalBuildFile(unittest.TestCase):


  def test_00(self):
    self.assertEqual( isGlobalBuildFile( 'Trilinos_version.h' ), True )


  def test_01(self):
    self.assertEqual( isGlobalBuildFile( 'CMakeLists.txt' ), True )


  def test_02(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/TrilinosPackages.cmake' ), True )


  def test_03(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/TrilinosCMakeQuickstart.txt' ), False )


  def test_04(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/ctest/experimental_build_test.cmake' ),
      False )


  def test_05(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/DependencyUnitTests/blah' ),
      False )


  def test_06(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/TPLs/FindTPLBLAS.cmake' ),
      True )


  def test_07(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/TPLs/FindTPLLAPACK.cmake' ),
      True )


  def test_08(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/TPLs/FindTPLMPI.cmake' ),
      True )


  def test_09(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/TPLs/FindTPLDummy.cmake' ),
      False )


  def test_10(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/utils/SetNotFound.cmake' ),
      True )


#
# Test extractPackageEnablesFromChangeStatus
#

class test_extractPackageEnablesFromChangeStatus(unittest.TestCase):


  def test_1(self):

    updateOutputStr = """
? packages/tpetra/doc/html
? packages/trilinoscouplings/doc/html
? packages/triutils/doc/html
? sampleScripts/checkin-test-gabriel.sh
M	cmake/TrilinosPackages.cmake
M	cmake/python/checkin-test.py
M	cmake/python/dump-cdash-deps-xml-file.py
P packages/thyra/dummy.blah
A	packages/teuchos/example/ExplicitInstantiation/four_files/CMakeLists.txt
"""

    options = MockOptions()
    enablePackagesList = []

    extractPackageEnablesFromChangeStatus(updateOutputStr, options,
      enablePackagesList, False)

    self.assertEqual( options.enableAllPackages, 'on' )
    self.assertEqual( enablePackagesList, ['Teuchos'] )


  def test_2(self):

    updateOutputStr = """
? packages/triutils/doc/html
M	cmake/python/checkin-test.py
M	cmake/python/dump-cdash-deps-xml-file.py
A	packages/nox/src/dummy.C
P packages/stratimikos/dummy.blah
M	packages/thyra/src/Thyra_ConfigDefs.hpp
M	packages/thyra/CMakeLists.txt
"""

    options = MockOptions()
    enablePackagesList = []

    extractPackageEnablesFromChangeStatus(updateOutputStr, options,
      enablePackagesList, False)

    self.assertEqual( options.enableAllPackages, 'default' )
    self.assertEqual( enablePackagesList, ['NOX', 'Thyra'] )

# Unit test class


class test_getLastCommitMessageStrFromRawCommitLogStr(unittest.TestCase):


  def test_01(self):
    cleanCommitMsg_expected = \
"""Some Commit Message

Some commit body

Some other message
"""
    rawLogOutput = "Standard git header stuff\n\n"+cleanCommitMsg_expected
    (cleanCommitMsg, numBlankLines) = getLastCommitMessageStrFromRawCommitLogStr(rawLogOutput)
    self.assertEqual(numBlankLines, -1)
    self.assertEqual(cleanCommitMsg, cleanCommitMsg_expected)


  def test_02(self):
    cleanCommitMsg_expected = \
"""Some Commit Message

Some commit body

Some other message
"""
    rawLogOutput = \
       "Standard git header stuff\n\n" \
       +cleanCommitMsg_expected+ \
       "\n\n\n=====================\n" \
       "Build/Test Cases Summary\n" \
       "=====================\n"
    (cleanCommitMsg, numBlankLines) = getLastCommitMessageStrFromRawCommitLogStr(rawLogOutput)
    self.assertEqual(numBlankLines, 3)
    self.assertEqual(cleanCommitMsg, cleanCommitMsg_expected)


  def test_03(self):
    cleanCommitMsg_expected = \
"""Some Commit Message

Some commit body

Some other message
"""
    rawLogOutput = \
       "Standard git header stuff\n\n" \
       +cleanCommitMsg_expected+ \
       "\n=====================\n" \
       "Build/Test Cases Summary\n" \
       "=====================\n"
    (cleanCommitMsg, numBlankLines) = getLastCommitMessageStrFromRawCommitLogStr(rawLogOutput)
    self.assertEqual(numBlankLines, 1)
    self.assertEqual(cleanCommitMsg, cleanCommitMsg_expected)


  def test_04(self):
    cleanCommitMsg_expected = \
"""Some Commit Message

Some commit body

Some other message"""
    rawLogOutput = \
       "Standard git header stuff\n\n" \
       +cleanCommitMsg_expected+ \
       "\n=====================\n" \
       "Build/Test Cases Summary\n" \
       "=====================\n"
    self.assertRaises(Exception, getLastCommitMessageStrFromRawCommitLogStr, rawLogOutput)


  def test_05(self):
    cleanCommitMsg_expected = \
"""Some Commit Message

Some commit body

Some other message
"""
    rawLogOutput = \
       "Standard git header stuff\n\n" \
       +cleanCommitMsg_expected+ \
       "\n=====================\n" \
       "Build/Test Cases Summary\n" \
       "=====================\n"
    (cleanCommitMsg, numBlankLines) = getLastCommitMessageStrFromRawCommitLogStr(rawLogOutput)
    self.assertEqual(numBlankLines, 1)
    self.assertEqual(cleanCommitMsg, cleanCommitMsg_expected)
    # Strip it again to make sure we can pull it off again and recover
    rawLogOutput = \
       "Standard git header stuff\n\n" \
       +cleanCommitMsg+ \
       "\n=====================\n" \
       "Build/Test Cases Summary\n" \
       "=====================\n"
    (cleanCommitMsg, numBlankLines) = getLastCommitMessageStrFromRawCommitLogStr(rawLogOutput)
    self.assertEqual(numBlankLines, 1)
    self.assertEqual(cleanCommitMsg, cleanCommitMsg_expected)


#
# Test getLastCommitMessageStrFromRawCommitLogStr
#

g_verbose=True
g_verbose=False

g_checkin_test_tests_dir = "checkin_test_tests"


def create_checkin_test_case_dir(testName, verbose=False):
  baseDir = os.getcwd()
  testDirName = os.path.join(g_checkin_test_tests_dir, testName)
  createDir(g_checkin_test_tests_dir, verbose)
  createDir(testDirName, verbose)
  return testDirName


# Unit test driver


def checkin_test_run_case(testObject, testName, optionsStr, cmndInterceptsStr, \
  expectPass, passRegexStrList, fromScratch=True, mustHaveCheckinTestOut=True \
  ):

  scriptsDir = getScriptBaseDir()
  verbose = g_verbose

  passRegexList = passRegexStrList.split('\n')

  if verbose: print "\npassRegexList =", passRegexList

  # A) Create the test directory

  baseDir = os.getcwd()
  echoChDir(create_checkin_test_case_dir(testName, verbose), verbose)

  try:

    # B) Create the command to run the checkin-test.py script
    
    cmnd = scriptsDir + "/checkin-test.py --no-eg-git-version-check " + optionsStr
    # NOTE: Above, we want to turn off the eg/git version tests since we want
    # these unit tests to run on machines that do not have the official
    # versions (e.g. the SCICO LAN) but where the versions might be okay.
    
    # C) Set up the command intercept file

    baseCmndInterceptsStr = \
      "FT: .*checkin-test-impl\.py.*\n" \
      "FT: which eg\n" \
      "FT: eg --version\n" \
      "FT: date\n" \
      "FT: rm [a-zA-Z0-9_/\.]+\n" \
      "FT: touch .*\n" \
      "FT: chmod .*\n" \
      "FT: hostname\n" \
      "FT: grep .* "+getTestOutputFileName()+"\n" \
      "FT: grep .*REQUESTED ACTIONS. PASSED.*\n"

    fullCmndInterceptsStr = baseCmndInterceptsStr + cmndInterceptsStr

    fullCmndInterceptsFileName = os.path.join(os.getcwd(), "cmndIntercepts.txt")
    writeStrToFile(fullCmndInterceptsStr, fullCmndInterceptsFileName)

    os.environ['GENERAL_SCRIPT_SUPPORT_CMND_INTERCEPTS_FILE'] = fullCmndInterceptsFileName
    
    # D) Run the checkin-test.py script with mock commands

    checkin_test_test_out = "checkin-test.test.out"

    rtnCode = echoRunSysCmnd(cmnd, timeCmnd=True, throwExcept=False,
      outFile=checkin_test_test_out, verbose=verbose)
    
    # E) Grep the output looking for specific string

    if mustHaveCheckinTestOut:
      outputFileToGrep = "checkin-test.out"
    else:
      outputFileToGrep = checkin_test_test_out
      
    for passRegex in passRegexList:
      foundRegex = getCmndOutput("grep '"+passRegex+"' "+outputFileToGrep, True, False)
      if verbose or not foundRegex:
        print "\ncheckin_test::"+testName+": Look for regex '"+passRegex+"' ...", 
        print "'"+foundRegex+"'", 
        if foundRegex: print ": PASSED"
        else: print ": FAILED"
      testObject.assertNotEqual(foundRegex, "")

    # F) Examine the final return code

    if expectPass:
      testObject.assertEqual(rtnCode, 0)
    else:
      testObject.assertNotEqual(rtnCode, 0)
    
  finally:
    # F) Get back to the current directory and reset
    echoChDir(baseDir, verbose=verbose)
    os.environ['GENERAL_SCRIPT_SUPPORT_CMND_INTERCEPTS_FILE']=""


# Test Data


g_cmndinterceptsInitialCommitPasses = \
  "IT: eg commit -a -F .*; 0; 'initial eg commit passes'\n"

g_cmndinterceptsDiffOnlyPasses = \
  "IT: eg diff --name-status.*; 0; 'M\tpackages/teuchos/CMakeLists.txt'\n"

g_cmndinterceptsPullOnlyPasses = \
  "IT: eg status; 1; 'eg status shows no uncommitted files'\n" \
  "IT: eg pull --rebase; 0; 'eg pull passed'\n"

g_cmndinterceptsPullPasses = \
  g_cmndinterceptsPullOnlyPasses \
  +g_cmndinterceptsDiffOnlyPasses

g_cmndinterceptsConfigPasses = \
  "IT: \./do-configure; 0; 'do-configure passed'\n"

g_cmndinterceptsConfigBuildPasses = \
  g_cmndinterceptsConfigPasses+ \
  "IT: make -j3; 0; 'make passed'\n"

g_cmndinterceptsConfigBuildTestPasses = \
  g_cmndinterceptsConfigBuildPasses+ \
  "IT: ctest -j5; 0; '100% tests passed, 0 tests failed out of 100'\n"

g_cmndinterceptsFinalPushPasses = \
  "IT: eg pull --rebase; 0; 'final eg pull --rebase passed'\n" \
  "IT: eg log --oneline origin..; 0; 'Only one commit'\n" \
  "IT: eg cat-file -p HEAD; 0; 'This is the last commit message'\n" \
  "IT: eg commit --amend -F .*; 0; 'Ammending the last commit'\n" \
  "IT: eg push; 0; 'push passes'\n"

g_cmndinterceptsSendBuildTestCaseEmail = \
  "IT: mailx -s .*; 0; 'Do not really sending build/test case email '\n"

g_cmndinterceptsSendFinalEmail = \
  "IT: sleep .*; 0; 'Do not really sleep'\n" \
  "IT: mailx -s .*; 0; 'Do not really send email '\n"

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
  "Test passed!\n" \
  "testResultsLine = 100% tests passed, 0 tests failed out of 100\n" \
  "passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
  "Test: Passed\n"

g_expectedRegexTestNotRun = \
  "The tests where never even run!\n" \
  "Test: FAILED\n"

g_expectedCommonOptionsSummary = \
  "Enabled Packages: Teuchos\n" \
  "Make Options: -j3\n" \
  "CTest Options: -j5\n"


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
    g_cmndinterceptsPullPasses \
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
    +g_expectedRegexTestPasses+ \
    "0) MPI_DEBUG => passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
    "1) SERIAL_RELEASE => Test case SERIAL_RELEASE was not run! => Does not affect commit/push readiness!\n" \
    +g_expectedCommonOptionsSummary+ \
    "=> A COMMIT IS OKAY TO BE PERFORMED!\n" \
    "=> A PUSH IS READY TO BE PERFORMED!\n" \
    "^READY TO PUSH: Trilinos:\n"
    )



#
# Test checkin_test
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
      "Usage: checkin-test.py \[OPTIONS\]\n" \
      "Quickstart:\n" \
      "Detailed Documentation:\n" \
      ".*--show-defaults.*\n" \
      ,
      mustHaveCheckinTestOut=False
      )
    # Help should not write the checkin-test.out file!
    self.assertEqual(
      os.path.exists(create_checkin_test_case_dir(testName, g_verbose)+"/checkin-test.out"),
      False)


  def test_do_all_commit_push_pass(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_commit_push_pass",
      "--make-options=-j3 --ctest-options=-j5" \
      " --commit-msg-header-file=cmake/python/utils/checkin_message_dummy1" \
      " --do-all --commit --push",
      \
      g_cmndinterceptsInitialCommitPasses \
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
      g_expectedRegexUpdateWithBuildCasePasses+ \
      g_expectedRegexConfigPasses+ \
      g_expectedRegexBuildPasses+ \
      g_expectedRegexTestPasses+ \
      "0) MPI_DEBUG => passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
      "1) SERIAL_RELEASE => passed: Trilinos/SERIAL_RELEASE: passed=100,notpassed=0\n" \
      +g_expectedCommonOptionsSummary+ \
      "=> A PUSH IS READY TO BE PERFORMED!\n" \
      "^DID PUSH: Trilinos:\n"
      )


  def test_do_all_without_serial_release_pass(self):
    g_test_do_all_without_serial_release_pass(self, "do_all_without_serial_release_pass")


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
      g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsFinalPushPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "0) MPI_DEBUG => passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
      "=> A PUSH IS READY TO BE PERFORMED!\n" \
      "^DID PUSH: Trilinos:\n"
      )


  def test_do_all_without_serial_release_then_commit_push_pass(self):

    testName = "do_all_without_serial_release_then_commit_push_pass"

    # Do the build/test only first (ready to push)
    g_test_do_all_without_serial_release_pass(self, testName)

    # Do the push after the fact
    checkin_test_run_case(
      \
      self,
      \
      testName,
      \
      "--make-options=-j3 --ctest-options=-j5 --without-serial-release --commit --push" \
      " --commit-msg-header-file=cmake/python/utils/checkin_message_dummy1" \
      ,
      \
      g_cmndinterceptsInitialCommitPasses \
      +g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsFinalPushPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "0) MPI_DEBUG => passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
      "=> A PUSH IS READY TO BE PERFORMED!\n" \
      "^DID PUSH: Trilinos:\n"
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
      "IT: eg diff --name-status.*; 0; 'eg diff passed'\n" 
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "0) MPI_DEBUG => passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
      "=> A PUSH IS READY TO BE PERFORMED!\n" \
      "^READY TO PUSH: Trilinos:\n"
      )


  def test_local_do_all_without_serial_release_pass(self):
    checkin_test_run_case(
      \
      self,
      \
      "local_do_all_without_serial_release_pass",
      \
      "--make-options=-j3 --ctest-options=-j5 --without-serial-release --local-do-all",
      \
      g_cmndinterceptsDiffOnlyPasses \
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
      +"0) MPI_DEBUG => passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
      +"1) SERIAL_RELEASE => Test case SERIAL_RELEASE was not run! => Does not affect commit/push readiness!\n" \
      +g_expectedCommonOptionsSummary \
      +"=> A COMMIT IS OKAY TO BE PERFORMED!\n" \
      +"A current successful pull does \*not\* exist => Not ready for final push!\n" \
      +"Explanation: In order to safely push, the local working directory needs to be up-to-date\n" \
      +"A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      +"^NOT READY TO PUSH: Trilinos:\n"
      )


  def test_do_all_without_serial_release_test_fail_force_commit_push_pass(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_without_serial_release_test_fail_force_commit_push_pass",
      \
      "--make-options=-j3 --ctest-options=-j5 --without-serial-release" \
      " --commit-msg-header-file=cmake/python/utils/checkin_message_dummy1" \
      " --do-all --commit --force-commit --push",
      \
      g_cmndinterceptsInitialCommitPasses \
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
      +"testResultsLine = 80% tests passed, 20 tests failed out of 100\n" \
      +"0) MPI_DEBUG => FAILED: Trilinos/MPI_DEBUG: passed=80,notpassed=20\n" \
      +"1) SERIAL_RELEASE => Test case SERIAL_RELEASE was not run! => Does not affect commit/push readiness!\n" \
      +g_expectedCommonOptionsSummary \
      +"Test: FAILED\n" \
      +"=> A COMMIT IS \*NOT\* OKAY TO BE PERFORMED!\n" \
      +"=> A PUSH IS READY TO BE PERFORMED!\n" \
      +"\*\*\* WARNING: The acceptance criteria for doing a commit/push has \*not\*\n" \
      +"\*\*\* been met, but a commit/push is being forced anyway by --force-commit!\n" \
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
      "FT: rm -rf MPI_DEBUG\n" \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "Running: rm -rf MPI_DEBUG\n" \
      +"0) MPI_DEBUG => No configure, build, or test for MPI_DEBUG was requested! => Not ready for final commit/push!\n" \
      +"=> A COMMIT IS \*NOT\* OKAY TO BE PERFORMED!\n" \
      +"=> A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      +"^NOT READY TO PUSH: Trilinos:\n"
      )


  # ToDo: Test --no-append-test-results
  # ToDo: Test --show-defaults
  
  
  # B) Test package enable/disable logic


  # ToDo: Test pulling info correctly from *.config files
  # ToDo: Test setting --enable-packages
  # ToDo: Test setting --disable-packages
  # ToDo: Test setting --enable-packages and --disable-packages
  # ToDo: Test setting --no-enable-fwd-packages
  # ToDo: Test setting --enable-all-packages=auto (changed from 'default')
  # ToDo: Test setting --enable-all-packages=off
  # ToDo: Test setting --enable-all-packages=on


  # C) Test partial actions short of running tests


  def test_without_serial_release_pull_only(self):
    checkin_test_run_case(
      self,
      \
      "without_serial_release_pull_only",
      \
      "--without-serial-release --pull",
      \
      g_cmndinterceptsPullPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      g_expectedRegexUpdatePasses \
      +"Not performing any build cases because no --configure, --build or --test was specified!\n" \
      +"0) MPI_DEBUG => No configure, build, or test for MPI_DEBUG was requested! => Not ready for final commit/push!\n" \
      +"A COMMIT IS \*NOT\* OKAY TO BE PERFORMED!\n" \
      +"A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      +"^NOT READY TO PUSH: Trilinos:\n"
      )


  def test_without_serial_release_pull_skip_commit_readiness_check(self):
    checkin_test_run_case(
      self,
      \
      "without_serial_release_pull_only",
      \
      "--without-serial-release --pull --skip-commit-readiness-check",
      \
      g_cmndinterceptsPullPasses \
      ,
      \
      True,
      \
      g_expectedRegexUpdatePasses \
      +"Skipping commit readiness check on request!\n" \
      +"Not performing commit/push or sending out commit/push readiness status on request!\n" \
      "^NOT READY TO PUSH$\n" \
      +"REQUESTED ACTIONS: PASSED\n"
      )


  def test_without_serial_release_pull_extra_pull_only(self):
    checkin_test_run_case(
      self,
      \
      "without_serial_release_pull_extra_pull_only",
      \
      "--pull --extra-pull-from='machine:/repo/dir/repo master'",
      \
      g_cmndinterceptsPullOnlyPasses \
      +"IT: eg pull --rebase machine:/repo/dir/repo master; 0; 'eg extra pull passed'\n"
      +g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      g_expectedRegexUpdatePasses \
      +"Pulling in updates from .machine:\/repo\/dir\/repo master.\n" \
      +"Not performing any build cases because no --configure, --build or --test was specified!\n" \
      +"A COMMIT IS \*NOT\* OKAY TO BE PERFORMED!\n" \
      +"A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      +"^NOT READY TO PUSH: Trilinos:\n"
      )


  def test_without_serial_release_configure_only(self):
    checkin_test_run_case(
      self,
      \
      "without_serial_release_configure_only",
      \
      "--without-serial-release --pull --configure",
      \
      g_cmndinterceptsPullPasses \
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
      "0) MPI_DEBUG => passed: Trilinos/MPI_DEBUG: configure-only passed => Not ready for final commit/push!\n" \
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
      g_cmndinterceptsPullPasses \
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
      "0) MPI_DEBUG => passed: Trilinos/MPI_DEBUG: build-only passed => Not ready for final commit/push!\n" \
      "A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      "NOT READY TO PUSH: Trilinos:\n"
      )


  # D) Test intermediate states with rerunning to fill out


  # ToDo: On all of these check that the right files are being deleted!


  # ToDo: Add test for pull followed by configure
  # ToDo: Add test for configure followed by build
  # ToDo: Add test for build followed by test
  # ToDo: Add test for removing files on pull
  # ToDo: Add test for removing files on configure
  # ToDo: Add test for removing files on build
  # ToDo: Add test for removing files on test


  # E) Test various failing use cases


  def test_do_all_no_eg_installed(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_no_eg_installed",
      \
      "--do-all" \
      ,
      \
      "IT: which eg; 1; '/usr/bin/which: no eg in (path1:path2:path3)'\n" \
      ,
      \
      False,
      \
      "Error, the eg command is not in your path!\n" \
      )


  def test_do_all_wrong_eg_vesion(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_wrong_eg_vesion",
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
      )


  def test_wrong_eg_vesion_ignore(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_wrong_eg_vesion",
      \
      "--no-eg-git-version-check --skip-commit-readiness-check" \
      ,
      \
      "IT: eg --version; 1; 'eg version wrong-version'\n" \
      ,
      \
      True,
      \
      "WARNING: No actions were performed!\n" \
      "REQUESTED ACTIONS: PASSED\n" \
      )

  #NOTE: I would also like to check the git verion but I can't becuase my
  #command intercept system can't hanlde more than one line of output.


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


  def test_do_all_commit_no_commit_msg_header(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_commit_no_commit_msg_header",
      \
      "--do-all --commit" \
      ,
      \
      "" \
      ,
      \
      False,
      \
      "Error, if you use --commit you must also specify --commit-msg-header-file!\n" \
      )


  def test_do_all_without_serial_release_commit_initial_commit_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_without_serial_release_commit_initial_commit_fail",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      " --commit-msg-header-file=cmake/python/utils/checkin_message_dummy1" \
      " --do-all --without-serial-release --commit" \
      ,
      \
      "IT: eg commit -a -F .*; 1; 'initial commit failed'\n" \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      "FAILED: Commit failed!\n" \
      "Commit failed, aborting pull!\n" \
      "Skipping getting list of modified files because pull failed!\n" \
      "The commit failed, skipping running the build/test cases!\n" \
      "0) MPI_DEBUG => The directory MPI_DEBUG does not exist! => Not ready for final commit/push!\n" \
      "A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      "Not doing the push on request (--no-push) but sending an email about the commit/push readiness status ...\n" \
      "INITIAL COMMIT FAILED: Trilinos:\n"
      )


  def test_do_all_without_serial_release_pull_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_without_serial_release_pull_fail",
      \
      "--do-all",
      "IT: eg status; 1; 'eg status shows no uncommitted files'\n" \
      "IT: eg pull --rebase; 1; 'eg pull failed'\n" \
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


  def test_do_all_without_serial_release_configure_fail(self):
    checkin_test_run_case(
      self,
      \
      "do_all_without_serial_release_configure_fail",
      "--do-all --without-serial-release",
      \
      g_cmndinterceptsPullPasses+ \
      "IT: \./do-configure; 1; 'do-configure failed'\n" \
      ,
      \
      False,
      \
      g_expectedRegexUpdateWithBuildCasePasses+ \
      "Configure failed returning 1!\n" \
      "The configure FAILED!\n" \
      "The build was never attempted!\n" \
      "The tests where never even run!\n" \
      "FAILED: Trilinos/MPI_DEBUG: configure failed\n" \
      "A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      "NOT READY TO PUSH: Trilinos:\n"
      )


  def test_do_all_without_serial_release_build_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_without_serial_release_build_fail",
      "--do-all --without-serial-release --make-options=-j3 --ctest-options=-j5",
      \
      g_cmndinterceptsPullPasses+ \
      "IT: \./do-configure; 0; 'do-configure passed'\n" \
      "IT: make -j3; 1; 'make filed'\n" \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      g_expectedRegexUpdateWithBuildCasePasses+ \
      g_expectedRegexConfigPasses+ \
      g_expectedRegexTestNotRun+ \
      g_expectedRegexBuildFailed+ \
      "0) MPI_DEBUG => FAILED: Trilinos/MPI_DEBUG: build failed => Not ready for final commit/push!\n" \
      "1) SERIAL_RELEASE => Test case SERIAL_RELEASE was not run! => Does not affect commit/push readiness!\n" \
      +g_expectedCommonOptionsSummary+ \
      "A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      "NOT READY TO PUSH: Trilinos:\n"
      )


  def test_do_all_without_serial_release_test_fail(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_without_serial_release_test_fail",
      "--do-all --without-serial-release --make-options=-j3 --ctest-options=-j5",
      \
      g_cmndinterceptsPullPasses+ \
      "IT: \./do-configure; 0; 'do-configure passed'\n" \
      "IT: make -j3; 0; 'make passed'\n" \
      "IT: ctest -j5; 1; '80% tests passed, 20 tests failed out of 100'\n" \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsSendFinalEmail \
      ,      \
      False,
      \
      g_expectedRegexUpdateWithBuildCasePasses+ \
      g_expectedRegexConfigPasses+ \
      g_expectedRegexBuildPasses+ \
      "FAILED: ctest failed returning 1!\n" \
      "testResultsLine = 80% tests passed, 20 tests failed out of 100\n" \
      "0) MPI_DEBUG => FAILED: Trilinos/MPI_DEBUG: passed=80,notpassed=20\n" \
      "1) SERIAL_RELEASE => Test case SERIAL_RELEASE was not run! => Does not affect commit/push readiness!\n" \
      +g_expectedCommonOptionsSummary+ \
      "Test: FAILED\n" \
      "A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      "NOT READY TO PUSH: Trilinos:\n"
      )


  # ToDo: Add test for the final pull failing
  # ToDo: Add test for final commit failing
  # ToDo: Add test for final push failing
  # ToDo: Add test that fails to push if some active pull was not performed
  #       first (in case --local-do-all is used).


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(testCheckinTest))
    return suite


if __name__ == '__main__':
  if os.path.exists(g_checkin_test_tests_dir):
    echoRunSysCmnd("rm -rf "+g_checkin_test_tests_dir, verbose=g_verbose)
  unittest.main()
