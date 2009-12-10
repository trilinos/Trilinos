########################################
# Unit testing code for CheckinTest.py #
########################################


from CheckinTest import *
import unittest


class MockOptions:
  def __init__(self):
    self.enableAllPackages = 'auto'


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
M	CMakeLists.txt
M	cmake/TrilinosPackages.cmake
M	cmake/python/checkin-test.py
M	 doc/Thyra/coding_guildlines/ThyraCodingGuideLines.tex
P packages/thyra/dummy.blah
A	packages/teuchos/example/ExplicitInstantiation/four_files/CMakeLists.txt
"""

    options = MockOptions()
    enablePackagesList = []

    extractPackageEnablesFromChangeStatus(updateOutputStr, options,
      enablePackagesList, False)

    self.assertEqual( options.enableAllPackages, 'on' )
    self.assertEqual( enablePackagesList, [u'TrilinosFramework', u'Teuchos'] )


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

    self.assertEqual( options.enableAllPackages, 'auto' )
    self.assertEqual( enablePackagesList, [u'TrilinosFramework', u'NOX', u'Thyra'] )



#############################################################################
#
#         Test getLastCommitMessageStrFromRawCommitLogStr
#
#############################################################################


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


#############################################################################
#
#                     Test checkin-test.py script
#
#############################################################################


# Test Data


g_cmndinterceptsInitialCommitPasses = \
  "IT: eg commit -a -F .*; 0; 'initial eg commit passed'\n"

g_cmndinterceptsDiffOnlyPasses = \
  "IT: eg diff --name-status.*; 0; 'M\tpackages/teuchos/CMakeLists.txt'\n"

g_cmndinterceptsPullOnlyPasses = \
  "IT: eg status; 1; 'eg status shows no uncommitted files'\n" \
  "IT: eg pull --rebase; 0; 'initial eg pull passed'\n"

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
  "IT: eg commit --amend -F .*; 0; 'Amending the last commit passed'\n" \
  "IT: eg push; 0; 'push passes'\n"

g_cmndinterceptsFinalPushNoAppendTestResultsPasses = \
  "IT: eg pull --rebase; 0; 'final eg pull --rebase passed'\n" \
  "IT: eg log --oneline origin..; 0; 'Only one commit'\n" \
  "IT: eg push; 0; 'push passes'\n"

g_cmndinterceptsSendBuildTestCaseEmail = \
  "IT: mailx -s .*; 0; 'Do not really sending build/test case email'\n"

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
  "No tests failed!\n" \
  "testResultsLine = .100% tests passed, 0 tests failed out of 100.\n" \
  "passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
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
  failRegexStrList=None, fileFailRegexStrList=None \
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
    
    cmnd = scriptsDir + "/../../checkin-test.py --no-eg-git-version-check " + optionsStr
    # NOTE: Above, we want to turn off the eg/git version tests since we want
    # these unit tests to run on machines that do not have the official
    # versions (e.g. the SCICO LAN) but where the versions might be okay.
    
    # C) Set up the command intercept file

    baseCmndInterceptsStr = \
      "FT: .*checkin-test-impl\.py.*\n" \
      "FT: eg config --get user.email\n" \
      "FT: which eg\n" \
      "FT: eg --version\n" \
      "FT: date\n" \
      "FT: rm [a-zA-Z0-9_/\.]+\n" \
      "FT: touch .*\n" \
      "FT: chmod .*\n" \
      "FT: hostname\n" \
      "FT: grep .* "+getTestOutputFileName()+"\n" \
      "FT: grep .*REQUESTED ACTIONS\: PASSED.*\n"

    fullCmndInterceptsStr = baseCmndInterceptsStr + cmndInterceptsStr

    fullCmndInterceptsFileName = os.path.join(os.getcwd(), "cmndIntercepts.txt")
    writeStrToFile(fullCmndInterceptsFileName, fullCmndInterceptsStr)

    os.environ['GENERAL_SCRIPT_SUPPORT_CMND_INTERCEPTS_FILE'] = fullCmndInterceptsFileName
    
    # D) Run the checkin-test.py script with mock commands

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
    +g_expectedRegexTestPasses \
    +"0) MPI_DEBUG: Will attempt to run!\n" \
    +"1) SERIAL_RELEASE: Will \*not\* attempt to run on request!\n" \
    +"0) MPI_DEBUG => passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
    +"1) SERIAL_RELEASE => Test case SERIAL_RELEASE was not run! => Does not affect commit/push readiness!\n" \
    +g_expectedCommonOptionsSummary \
    +"=> A COMMIT IS OKAY TO BE PERFORMED!\n" \
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
    " --allow-no-pull --configure --send-email-to= --skip-commit-readiness-check" \
    +" " +optionsStr \
    ,
    \
    "IT: eg diff --name-status.*; 0; '"+modifiedFilesStr+"'\n" \
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


  def test_do_all_commit_push_pass(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_commit_push_pass",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      +" --commit-msg-header-file=cmake/python/utils/checkin_message_dummy1" \
      +" --do-all --commit --push" \
      +" --execute-on-ready-to-push=\"ssh -q godel /some/dir/some_command.sh &\"",
      \
      g_cmndinterceptsInitialCommitPasses \
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
      g_expectedRegexUpdateWithBuildCasePasses \
      +g_expectedRegexConfigPasses \
      +g_expectedRegexBuildPasses \
      +g_expectedRegexTestPasses \
      +"0) MPI_DEBUG => passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
      +"1) SERIAL_RELEASE => passed: Trilinos/SERIAL_RELEASE: passed=100,notpassed=0\n" \
      +g_expectedCommonOptionsSummary \
      +"=> A PUSH IS READY TO BE PERFORMED!\n" \
      +"mailx .* trilinos-checkin-tests.*\n" \
      +"^DID PUSH: Trilinos:\n" \
      +"Executing final command (ssh -q godel /some/dir/some_command.sh &) since a push is okay to be performed!\n" \
      +"Running: ssh -q godel /some/dir/some_command.sh &\n" \
      ,
      [
      (getStatusOutputFileName(), "eg status shows no uncommitted files\n"),
      (getInitialCommitOutputFileName(), "initial eg commit passed\n"),
      (getInitialPullOutputFileName(), "initial eg pull passed\n"),
      (getModifiedFilesOutputFileName(), "M\tpackages/teuchos/CMakeLists.txt\n"),
      (getFinalPullOutputFileName(), "final eg pull --rebase passed\n"),
      (getFinalCommitEmailBodyFileName(),
         getAutomatedStatusSummaryHeaderKeyStr()+"\n"
         +"Enabled Packages: Teuchos\n" \
         +"Enabled all Forward Packages\n" \
         ),
      ("MPI_DEBUG/do-configure.base",
       "\-DTPL_ENABLE_MPI:BOOL=ON\n" \
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
       "\-DTrilinos_ENABLE_TESTS:BOOL=ON\n" \
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
      +"^NOT READY TO PUSH: Trilinos:\n" \
      +"Not executing final command (ssh -q godel /some/dir/some_command.sh &) since a push is not okay to be performed!\n" \
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
      " --do-all --commit --force-commit-push --push",
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
      +"testResultsLine = .80% tests passed, 20 tests failed out of 100.\n" \
      +"0) MPI_DEBUG => FAILED: Trilinos/MPI_DEBUG: passed=80,notpassed=20\n" \
      +"1) SERIAL_RELEASE => Test case SERIAL_RELEASE was not run! => Does not affect commit/push readiness!\n" \
      +g_expectedCommonOptionsSummary \
      +"Test: FAILED\n" \
      +"=> A COMMIT IS \*NOT\* OKAY TO BE PERFORMED!\n" \
      +"=> A PUSH IS READY TO BE PERFORMED!\n" \
      +"\*\*\* WARNING: The acceptance criteria for doing a commit/push has \*not\*\n" \
      +"\*\*\* been met, but a commit/push is being forced anyway by --force-commit-push!\n" \
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


  def test_do_all_no_append_test_results_push_pass(self):
    checkin_test_run_case(
      \
      self,
      \
      "do_all_no_append_test_results_push_pass",
      \
      "--make-options=-j3 --ctest-options=-j5" \
      +" --commit-msg-header-file=cmake/python/utils/checkin_message_dummy1" \
      +" --do-all --no-append-test-results --push",
      \
      g_cmndinterceptsPullPasses \
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
      +"0) MPI_DEBUG => passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
      +"1) SERIAL_RELEASE => passed: Trilinos/SERIAL_RELEASE: passed=100,notpassed=0\n" \
      +g_expectedCommonOptionsSummary \
      +"Skipping appending test results on request (--no-append-test-results)!\n" \
      +"=> A PUSH IS READY TO BE PERFORMED!\n" \
      +"^DID PUSH: Trilinos:\n" \
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
      "--enable-packages=TrilinosFramework,RTOp,Thyra" \
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
      )


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
      "--pull --extra-pull-from=machine:/repo/dir/repo:master",
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
      +"eg pull --rebase machine:\/repo\/dir\/repo master\n" \
      +"Not performing any build cases because no --configure, --build or --test was specified!\n" \
      +"A COMMIT IS \*NOT\* OKAY TO BE PERFORMED!\n" \
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
      g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      "Skipping all updates on request!\n" \
      +"Not performing any build cases because no --configure, --build or --test was specified!\n" \
      +"A COMMIT IS \*NOT\* OKAY TO BE PERFORMED!\n" \
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
      g_cmndinterceptsDiffOnlyPasses \
      +g_cmndinterceptsFinalPushPasses \
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      True,
      \
      "0) MPI_DEBUG => passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
      "0) MPI_DEBUG Results:\n" \
      +"=> A PUSH IS READY TO BE PERFORMED!\n" \
      +"^DID PUSH: Trilinos:\n" \
      
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
      +" --extra-pull-from=dummy:master" \
      +" --commit-msg-header-file=cmake/python/utils/checkin_message_dummy1" \
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
      +"=> A PUSH IS READY TO BE PERFORMED!\n" \
      +"^DID PUSH: Trilinos:\n" \
      ,
      failRegexStrList = \
      "eg pull --rebase dummy master\n" \
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


  # ToDo: On all of these below check that the right files are being deleted!

  # ToDo: Add test for removing files on pull (fail immediately)
  # ToDo: Add test for removing files on configure (fail immediately)
  # ToDo: Add test for removing files on build (fail immediately)
  # ToDo: Add test for removing files on test (fail immediately)


  # F) Test various failing use cases


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
      g_cmndinterceptsDiffOnlyPasses \
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
      +"0) MPI_DEBUG => FAILED: Trilinos/MPI_DEBUG: pre-configure failed => Not ready for final commit/push!\n" \
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
      g_cmndinterceptsPullPasses+ \
      "IT: \./do-configure; 1; 'do-configure failed'\n" \
      ,
      \
      False,
      \
      g_expectedRegexUpdateWithBuildCasePasses \
      +"Configure failed returning 1!\n" \
      +"The configure FAILED!\n" \
      +"The build was never attempted!\n" \
      +"The tests where never even run!\n" \
      +"FAILED: Trilinos/MPI_DEBUG: configure failed\n" \
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
      g_cmndinterceptsPullPasses+ \
      "IT: \./do-configure; 0; 'do-configure passed'\n" \
      "IT: make -j3; 1; 'make filed'\n" \
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
      +"0) MPI_DEBUG => FAILED: Trilinos/MPI_DEBUG: build failed => Not ready for final commit/push!\n" \
      +"1) SERIAL_RELEASE => Test case SERIAL_RELEASE was not run! => Does not affect commit/push readiness!\n" \
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
      g_cmndinterceptsPullPasses \
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
      +"0) MPI_DEBUG => FAILED: Trilinos/MPI_DEBUG: passed=80,notpassed=20\n" \
      +"1) SERIAL_RELEASE => Test case SERIAL_RELEASE was not run! => Does not affect commit/push readiness!\n" \
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
      g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +"IT: eg pull --rebase; 1; 'final eg pull --rebase FAILED'\n" \
      +"IT: eg log --oneline origin..; 0; 'Only one commit'\n" \
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
      g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +"IT: eg pull --rebase; 0; 'final eg pull --rebase passed'\n" \
      +"IT: eg log --oneline origin..; 0; 'Only one commit'\n" \
      +"IT: eg cat-file -p HEAD; 0; 'This is the last commit message'\n" \
      +"IT: eg commit --amend -F .*; 1; 'Amending the last commit FAILED'\n" \
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
      +"Attempting to ammend the final commmit message ...\n" \
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
      g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +"IT: eg pull --rebase; 0; 'final eg pull --rebase passed'\n" \
      +"IT: eg log --oneline origin..; 0; 'Only one commit'\n" \
      +"IT: eg cat-file -p HEAD; 0; 'This is the last commit message'\n" \
      +"IT: eg commit --amend -F .*; 0; 'Amending the last commit passed'\n" \
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
      g_cmndinterceptsPullPasses \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +g_cmndinterceptsConfigBuildTestPasses \
      +g_cmndinterceptsSendBuildTestCaseEmail \
      +"IT: eg pull --rebase; 0; 'final eg pull --rebase passed'\n" \
      +"IT: eg log --oneline origin..; 0; ''\n" \
      +"IT: eg cat-file -p HEAD; 0; 'Some commit not the local commit'\n" \
      +"IT: eg push; 1; 'push FAILED due to no local commits'\n"
      +g_cmndinterceptsSendFinalEmail \
      ,
      \
      False,
      \
      g_expectedRegexUpdateWithBuildCasePasses \
      +g_expectedRegexConfigPasses \
      +g_expectedRegexBuildPasses \
      +g_expectedRegexTestPasses \
      +"0) MPI_DEBUG => passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
      +"1) SERIAL_RELEASE => passed: Trilinos/SERIAL_RELEASE: passed=100,notpassed=0\n" \
      +"=> A PUSH IS READY TO BE PERFORMED!\n" \
      +"No local commits exit!\n" \
      +"Skipping ammending last commit because there are no local commits!\n" \
      +"Attempting to do the push ...\n" \
      +"Push failed!\n" \
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
      g_cmndinterceptsPullPasses \
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
       +"0) MPI_DEBUG => FAILED: Trilinos/MPI_DEBUG: no tests run\n" \
      +"=> A COMMIT IS \*NOT\* OKAY TO BE PERFORMED!\n" \
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
      g_cmndinterceptsDiffOnlyPasses \
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
      +"0) MPI_DEBUG => passed: Trilinos/MPI_DEBUG: passed=100,notpassed=0\n" \
      +"=> A COMMIT IS OKAY TO BE PERFORMED!\n" \
      +"A current successful pull does \*not\* exist => Not ready for final push!\n"\
      +"=> A PUSH IS \*NOT\* READY TO BE PERFORMED!\n"\
      +"^ABORTED COMMIT/PUSH: Trilinos:\n" \
      +"REQUESTED ACTIONS: FAILED\n" \
      )


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(testCheckinTest))
    return suite


if __name__ == '__main__':
  if os.path.exists(g_checkin_test_tests_dir):
    echoRunSysCmnd("rm -rf "+g_checkin_test_tests_dir, verbose=g_verbose)
  unittest.main()
