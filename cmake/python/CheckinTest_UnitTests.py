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


# Unit test driver


def checkin_test_run_case(testObject, testName, optionsStr, cmndInterceptsStr, \
  expectPass, passRegexStrList, fromScratch=True \
  ):

  scriptsDir = getScriptBaseDir()
  verbose = g_verbose

  passRegexList = passRegexStrList.split('\n')

  if verbose: print "\npassRegexList =", passRegexList

  # A) Create the test directory

  baseDir = os.getcwd()
  createDir(g_checkin_test_tests_dir, True, verbose)
  createDir(testName, True, verbose)

  try:

    # B) Create the command to run the checkin-test.py script
    
    cmnd = scriptsDir + "/checkin-test.py " + optionsStr
    
    # C) Set up the command intercept file

    baseCmndInterceptsStr = \
      "FT: .*checkin-test-impl\.py.*\n" \
      "FT: date\n" \
      "FT: rm .*\n" \
      "FT: touch .*\n" \
      "FT: chmod .*\n" \
      "FT: hostname\n" \
      "FT: grep .*OVERALL. PASSED.*\n"

    fullCmndInterceptsStr = baseCmndInterceptsStr + cmndInterceptsStr

    fullCmndInterceptsFileName = os.path.join(os.getcwd(), "cmndIntercepts.txt")
    writeStrToFile(fullCmndInterceptsStr, fullCmndInterceptsFileName)

    os.environ['GENERAL_SCRIPT_SUPPORT_CMND_INTERCEPTS_FILE'] = fullCmndInterceptsFileName
    
    # D) Run the checkin-test.py script with mock commands

    rtnCode = echoRunSysCmnd(cmnd, timeCmnd=True, throwExcept=False,
      outFile="checkin-test.test.out", verbose=verbose)
    
    # E) Grep the output looking for specific string

    for passRegex in passRegexList:
      foundRegex = getCmndOutput("grep '"+passRegex+"' checkin-test.out", True, False)
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


# ToDo: Extract out expected regex strings as helper varaibles and switch from
# an array to a single string and then split on '\n'.


g_cmndinterceptsPullPasses = \
  "IT: eg status; 1; 'eg status shows no uncommitted files'\n" \
  "IT: eg pull --rebase; 0; 'eg pull passed'\n" \
  "IT: eg diff --name-status.*; 0; 'M\tpackages/teuchos/CMakeLists.txt'\n"

g_cmndinterceptsConfigBuildTestPasses = \
  "IT: \./do-configure; 0; 'do-configure passed'\n" \
  "IT: make -j3; 0; 'make passed'\n" \
  "IT: ctest -j5; 0; '100% tests passed, 0 tests failed out of 100'\n"

g_cmndinterceptsFinalPullCommitPasses = \
  "IT: eg pull --rebase; 0; 'final eg pull --rebase passed'\n" \
  "IT: eg log --oneline origin..; 0; 'Only one commit'\n" \
  "IT: eg cat-file -p HEAD; 0; 'This is the last commit message'\n" \
  "IT: eg commit --amend -F .*; 0; 'Ammending the last commit'\n" \

g_cmndinterceptsSendEmail = \
  "IT: mailx -s .*; 0; 'Do not really send email '\n"

g_cmndinterceptsSendFinalEmail = \
  "IT: sleep .*; 0; 'Do not really sleep'\n" \
  "IT: mailx -s .*; 0; 'Do not really send email '\n"

g_expectedRegexUpdatePasses = \
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

#
# Test checkin_test
#

class test_checkin_test(unittest.TestCase):


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
      "FT: grep .*\n" \
      "IT: eg commit -a -F .*; 0; 'initial eg commit passes'\n" \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendEmail \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendEmail \
      +g_cmndinterceptsFinalPullCommitPasses+ \
      "IT: eg push; 0; 'push passes'\n" \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      g_expectedRegexUpdatePasses+ \
      g_expectedRegexConfigPasses+ \
      g_expectedRegexBuildPasses+ \
      g_expectedRegexTestPasses+ \
      "0) MPI_DEBUG => passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
      "1) SERIAL_RELEASE => passed: Trilinos/SERIAL_RELEASE: passed=100,notpassed=0\n" \
      +g_expectedCommonOptionsSummary+ \
      "=> A PUSH IS OKAY TO BE PERFORMED!\n" \
      "^DID PUSH: Trilinos:\n"
      )


  def test_do_all_without_serial_release_pass(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_without_serial_release_pass",
      "--do-all --without-serial-release --make-options=-j3 --ctest-options=-j5",
      \
      "FT: grep .*\n" \
      +g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendEmail \
      +g_cmndinterceptsFinalPullCommitPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      g_expectedRegexUpdatePasses \
      +g_expectedRegexConfigPasses \
      +g_expectedRegexBuildPasses \
      +g_expectedRegexTestPasses+ \
      "0) MPI_DEBUG => passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
      "1) SERIAL_RELEASE => Test case SERIAL_RELEASE was not run!  Does not affect commit/push readiness!\n" \
      +g_expectedCommonOptionsSummary+ \
      "=> A PUSH IS OKAY TO BE PERFORMED!\n" \
      "^READY TO PUSH: Trilinos:\n"
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
      g_expectedRegexUpdatePasses+ \
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
      +g_cmndinterceptsSendEmail \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      g_expectedRegexUpdatePasses+ \
      g_expectedRegexConfigPasses+ \
      g_expectedRegexTestNotRun+ \
      g_expectedRegexBuildFailed+ \
      "0) MPI_DEBUG => FAILED: Trilinos/MPI_DEBUG: build failed  Not ready for final commit/push!\n" \
      "1) SERIAL_RELEASE => Test case SERIAL_RELEASE was not run!  Does not affect commit/push readiness!\n" \
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
      "FT: grep .*\n" \
      +g_cmndinterceptsPullPasses+ \
      "IT: \./do-configure; 0; 'do-configure passed'\n" \
      "IT: make -j3; 0; 'make passed'\n" \
      "IT: ctest -j5; 1; '80% tests passed, 20 tests failed out of 100'\n" \
      +g_cmndinterceptsSendEmail \
      +g_cmndinterceptsSendFinalEmail \
      ,      \
      False,
      \
      g_expectedRegexUpdatePasses+ \
      g_expectedRegexConfigPasses+ \
      g_expectedRegexBuildPasses+ \
      "FAILED: ctest failed returning 1!\n" \
      "testResultsLine = 80% tests passed, 20 tests failed out of 100\n" \
      "0) MPI_DEBUG => FAILED: Trilinos/MPI_DEBUG: passed=80,notpassed=20\n" \
      "1) SERIAL_RELEASE => Test case SERIAL_RELEASE was not run!  Does not affect commit/push readiness!\n" \
      +g_expectedCommonOptionsSummary+ \
      "Test: FAILED\n" \
      "A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      "NOT READY TO PUSH: Trilinos:\n"
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
      "0) MPI_DEBUG => The directory MPI_DEBUG does not exist!  Not ready for final commit/push!\n" \
      "A PUSH IS \*NOT\* READY TO BE PERFORMED!\n" \
      "Not attempted final commit and/or push!\n" \
      "INITIAL COMMIT FAILED: Trilinos:\n"
      )


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(testCheckinTest))
    return suite


if __name__ == '__main__':
  if os.path.exists(g_checkin_test_tests_dir):
    echoRunSysCmnd("rm -rf "+g_checkin_test_tests_dir, verbose=g_verbose)
  unittest.main()
