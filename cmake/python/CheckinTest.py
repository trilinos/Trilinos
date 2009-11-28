
#
# ToDo:
#
#  (*) Create a TaskStatus class and use it to simplify the logic replacing
#  the simple bools.
#
#  (*) Put in checks for the names of Trilinos packages from --enable-packages
#  and --disable-packages arguments.  Right now a mispelled package name would
#  just be ignored.  Also, put in unit tests for this.
#
#  (*) Implement --extra-builds option.
#
#  (*) Implement check that -DTPL_ENABLE* is not specified in any of the
#  standard configure files.  Also, make sure that there are no
#  -DTrilinos_ENABLE* variables in COMMON.config either.  Force users to do
#  package enables/disables using --enable-packages, --disable-packages.
#
#  (*) Change logic to not enable everything if TrilinosPackages.cmake or
#  TrilinosTPLs.cmake are changed.
#
#  (*) Add --execute-on-pass to set a command to launch after everything in
#  the script has occured if the tests passed.  This can be used for launching
#  a remote test/push with ssh.  Remember to print out this command in the
#  final summary email to let user know that it is running.  Hint: Users can
#  run with:
#
#    --execute-on-pass="ssh -q godel 'checkin-test-godel.sh --do-all --extra-pull-from=... 2>&1 > /dev/null' &"
#
#  (*) Implement checking and control for enabling TrilinosFramework tests or
#  not depending on if anything under cmake/ or packages/teuchos/test/CTestxxx
#  is modified or not, except if TrilinosPackages.cmake and TrilinosTPLs.cmake
#  is specified.  Also, implement Teuchos in this case.
#
#  (*) Make this work automatically on branches too.  It should pull and push
#  to 'origin' always but the current branch should always be used.  Change
#  --extra-pull-from to require the repository and branch.
#
#  (*) Once everyone is using the checkin-test.py script:
#
#  - Turn off framework tests by default and turn them in checkin
#    testing ...
#
#  - Turn off generation of HTML/XML files by default and turn them on in
#    checkin testing ...
#

#
# General scripting support
#
# NOTE: Included first to check the version of python!
#

from GeneralScriptSupport import *

from TrilinosDependencies import getTrilinosDependenciesFromXmlFile
from TrilinosDependencies import defaultTrilinosDepsXmlInFile
from TrilinosPackageFilePathUtils import *
import time


# Set the official eg/git versions!
g_officialEgVersion = "1.6.5.3"
g_officialGitVersion = "1.6.5.2"


def getCommonConfigFileName():
  return "COMMON.config"


def getTestCaseName(serialOrMpi, buildType):
  return serialOrMpi + "_" + buildType

def getBuildSpecificConfigFileName(serialOrMpi, buildType):
  return getTestCaseName(serialOrMpi, buildType) + ".config"


def getStatusOutputFileName():
  return "status.out"


def getInitialPullOutputFileName():
  return "pullInitial.out"


def getInitialExtraPullOutputFileName():
  return "pullInitialExtra.out"


def getInitialPullSuccessFileName():
  return "pullInitial.success"


def getModifiedFilesOutputFileName():
  return "modifiedFiles.out"


def getFinalPullOutputFileName():
  return "pullFinal.out"


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


def getInitialCommitEmailBodyFileName():
  return "commitInitialEmailBody.out"


def getInitialCommitOutputFileName():
  return "commitInitial.out"


def getFinalCommitEmailBodyFileName():
  return "commitFinalEmailBody.out"


def getFinalCommitOutputFileName():
  return "commitFinal.out"


def getCommitStatusEmailBodyFileName():
  return "commitStatusEmailBody.out"


def getPushOutputFileName():
  return "push.out"


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
  

def assertEgGitVersions(inOptions):

  egWhich = getCmndOutput("which eg", True, False)
  if egWhich == "" or re.match(".+no eg.+", egWhich):
    raise Exception("Error, the eg command is not in your path! ("+egWhich+")")

  egVersionOuput = getCmndOutput("eg --version", True, False)
  egVersionsList = egVersionOuput.split('\n')

  if inOptions.enableEgGitVersionCheck:
    assertEgGitVersionHelper(egVersionsList[0], "eg version "+g_officialEgVersion)
    assertEgGitVersionHelper(egVersionsList[1], "git version "+g_officialGitVersion)


def executePull(inOptions, baseTestDir, outFile, pullFromRepo=None):
  cmnd = "eg pull --rebase"
  if pullFromRepo:
    print "\nPulling in updates from '"+pullFromRepo+"' ...\n"
    cmnd += " " + pullFromRepo
  else:
    print "\nPulling in updates from 'origin' ...\n"
  return echoRunSysCmnd( cmnd,
    workingDir=inOptions.trilinosSrcDir,
    outFile=os.path.join(baseTestDir, outFile),
    timeCmnd=True, returnTimeCmnd=True, throwExcept=False
    )
  

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


def writeDefaultBuildSpecificConfigFile(serialOrMpi, buildType):

  buildSpecificConfigFileName = getBuildSpecificConfigFileName(serialOrMpi, buildType)

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


def readAndAppendCMakeOptions(fileName, cmakeOptions_inout):

  if not os.path.exists(fileName):
    return

  print "\nAppending options from "+fileName+":"

  cmakeOptionsFile = open(fileName, 'r')

  for line in cmakeOptionsFile:
    if line[0] != '#':
      print "  Appnding: "+line.strip()
      cmakeOptions_inout.append(line.strip())


reModifedFiles = re.compile(r"^[MA]\t(.+)$")


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


def getCurrentDiffOutput(inOptions, baseTestDir):
  echoRunSysCmnd(
    "eg diff --name-status origin",
    workingDir=inOptions.trilinosSrcDir,
    outFile=os.path.join(baseTestDir, getModifiedFilesOutputFileName()),
    timeCmnd=True
    )


def extractPackageEnablesFromChangeStatus(updateOutputStr, inOptions_inout,
  enablePackagesList_inout, verbose=True ) \
  :

  trilinosDependencies = getTrilinosDependenciesFromXmlFile(defaultTrilinosDepsXmlInFile)

  modifiedFilesList = extractFilesListMatchingPattern(
    updateOutputStr.split('\n'), reModifedFiles )

  for modifiedFileFullPath in modifiedFilesList:

    if isGlobalBuildFile(modifiedFileFullPath):
      if inOptions_inout.enableAllPackages == 'default':
        if verbose:
          print "\nModifed file: '"+modifiedFileFullPath+"'\n" \
            "  => Enabling all Trilinos packages!"
        inOptions_inout.enableAllPackages = 'on'

    packageName = getPackageNameFromPath(trilinosDependencies, modifiedFileFullPath)
    if packageName and findInSequence(enablePackagesList_inout, packageName) == -1:
      if verbose:
        print "\nModified file: '"+modifiedFileFullPath+"'\n" \
          "  => Enabling '"+packageName+"'!"
      enablePackagesList_inout.append(packageName)


def createConfigureFile(cmakeOptions, baseCmnd, trilinosSrcDir, configFileName):

    doConfigStr = ""
  
    doConfigStr += \
      "EXTRA_ARGS=$@\n" \
      "\n" \
      +baseCmnd+ " \\\n"
  
    for opt in cmakeOptions:
      doConfigStr += opt + " \\\n"
    
    doConfigStr += \
      "$EXTRA_ARGS"

    if trilinosSrcDir:
      doConfigStr += " \\\n"+trilinosSrcDir
    
    doConfigStr += "\n"
  
    writeStrToFile(configFileName, doConfigStr)
    echoRunSysCmnd('chmod a+x '+configFileName)


reCtestFailTotal = re.compile(r".+, ([0-9]+) tests failed out of ([0-9]+)")


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


def getStageStatus(stageName, stageDoBool, stagePassed, stageTiming):
  stageStatusStr = stageName + ": "
  if stageDoBool:
    if stagePassed:
      stageStatusStr += "Passed"
    else:
      stageStatusStr += "FAILED"
    stageStatusStr += " ("+str(stageTiming)+" min)"
  else:
    stageStatusStr += "Not Performed"
  stageStatusStr += "\n"
  return stageStatusStr


def analyzeResultsSendEmail(inOptions, buildDirName,
  enabledPackagesList, cmakeOptions, startingTime, timings ) \
  :

  print ""
  print "1) Determine what passed and failed ..."
  print ""

  success = False

  # Determine if the update passed

  commitPassed = None
  updatePassed = None
  updateOutputExists = False

  if inOptions.doPull:

    if os.path.exists("../"+getInitialPullOutputFileName()):
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

    if os.path.exists(getTestSuccessFileName()):
      print "\nAll of the tests ran passed!\n"
      testsPassed = True
    elif testOutputExists:
      print "\nAt least one of the tests ran FAILED!\n"
      testsPassed = False
    else:
      print "\nThe tests where never even run!\n"
      testsPassed = False

    if testOutputExists:

      testResultsLine = getCmndOutput("grep 'tests failed out of' "+getTestOutputFileName(),
        True, False)

      print "testResultsLine =", testResultsLine

      reCtestFailTotalMatch = reCtestFailTotal.match(testResultsLine)

      if reCtestFailTotalMatch:
        numFailedTests = int(reCtestFailTotalMatch.group(1))
        numTotalTests = int(reCtestFailTotalMatch.group(2))
        numPassedTests = numTotalTests - numFailedTests
      else:
        numTotalTests = None
        numPassedTests = None
        testsPassed = False

  else:

    print "\nRunning the tests was not performed!\n"

  print ""
  print "2) Construct the email message ..."
  print ""

  # 2.a) Construct the subject line

  overallPassed = None
  subjectLine = "Trilinos/"+buildDirName
  selectedFinalStatus = False

  if inOptions.doTest and not selectedFinalStatus:
    if testOutputExists:
      if numTotalTests:
        subjectLine += ": passed="+str(numPassedTests)+",notpassed="+str(numFailedTests)
      else:
        subjectLine += ": no tests run"
      if testsPassed:
        overallPassed = True
      else:
        overallPassed = False
      selectedFinalStatus = True
    elif not inOptions.doBuild and not buildOutputExists:
      subjectLine += ": no active build exists"
      overallPassed = False
      selectedFinalStatus = True

  if inOptions.doBuild and not selectedFinalStatus:
    if buildPassed:
      subjectLine += ": build-only passed"
      overallPassed = True
      selectedFinalStatus = True
    elif buildOutputExists:
      subjectLine += ": build failed"
      overallPassed = False
      selectedFinalStatus = True

  if inOptions.doConfigure and not selectedFinalStatus:
    if configurePassed:
      subjectLine += ": configure-only passed"
      overallPassed = True
      selectedFinalStatus = True
    elif configureOutputExists:
      subjectLine += ": configure failed"
      overallPassed = False
      selectedFinalStatus = True
      selectedFinalStatus = True

  if inOptions.doPull and not selectedFinalStatus:
    if updatePassed:
      subjectLine += ": update-only passed"
      overallPassed = True
      selectedFinalStatus = True
    elif updateOutputExists:
      subjectLine += ": update FAILED"
      overallPassed = False
      selectedFinalStatus = True

  if not selectedFinalStatus:
    raise Exception("Error, final pass/fail status not found!")

  if overallPassed:
    subjectLine = "passed: " + subjectLine
  else:
    subjectLine = "FAILED: " + subjectLine

  print "\nsubjectLine = '"+subjectLine+"'\n"

  success = overallPassed

  # 2.b) Construct the email body

  emailBody = subjectLine + "\n\n"

  emailBody += getCmndOutput("date", True) + "\n\n"

  emailBody += "Enabled Packages: " + ', '.join(enabledPackagesList) + "\n"
  if inOptions.disablePackages:
    emailBody += "Disabled Packages: " + inOptions.disablePackages + "\n"
  emailBody += "Hostname: " + getHostname() + "\n"
  emailBody += "Source Dir: " + inOptions.trilinosSrcDir + "\n"
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

    if inOptions.showAllTests:
      emailBody += getCmndOutput("cat "+getTestOutputFileName())
    else:
      emailBody += getCmndOutput("grep -A 10000 '\% tests passed, ' "+getTestOutputFileName())

  else:

    emailBody += "\n***\n*** WARNING: There are no test results!\n***\n\n"

  endingTime = time.time()
  totalTime = (endingTime - startingTime) / 60.0

  emailBody += "\nTotal time for "+buildDirName+" = "+str(totalTime) + " minutes"

  #print "emailBody:\n\n\n\n", emailBody, "\n\n\n\n"

  writeStrToFile(getEmailBodyFileName(), emailBody)

  if overallPassed:
    echoRunSysCmnd("touch "+getEmailSuccessFileName())

  print ""
  print "3) Send the email message ..."
  print ""

  if inOptions.sendEmailTo:

    emailAddresses = getEmailAddressesSpaceString(inOptions.sendEmailTo)
    echoRunSysCmnd("mailx -s \""+subjectLine+"\" "+emailAddresses+" < "+getEmailBodyFileName())

  else:

    print "Not sending email because no email addresses were given!"

  # 3) Return final result

  return success


def getTestCaseSummaryLine(testCaseName):
  # Get the email file
  absEmailBodyFileName = testCaseName+"/"+getEmailBodyFileName()
  if os.path.exists(absEmailBodyFileName):
    testCaseEmailStrArray = open(absEmailBodyFileName, 'r').readlines()
  else:
    testCaseEmailStrArray = None
  # Get the first line (which is the summary)
  if testCaseEmailStrArray:
    summaryLine = testCaseEmailStrArray[0].strip()
  else:
    summaryLine = \
      "Error, The build/test was never completed!" \
      " (the file '"+absEmailBodyFileName+"' does not exist.)"
  return summaryLine


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


def getSummaryEmailSectionStr(inOptions):
  summaryEmailSectionStr = ""
  if performAnyBuildTestActions(inOptions):
    if inOptions.withMpiDebug:
      summaryEmailSectionStr += getTestCaseEmailSummary("MPI_DEBUG", 0)
    if inOptions.withSerialRelease:
      summaryEmailSectionStr += getTestCaseEmailSummary("SERIAL_RELEASE", 1)
  return summaryEmailSectionStr

  
def runTestCase(inOptions, serialOrMpi, buildType, extraCMakeOptions,
  timings ) \
  :

  success = True

  startingTime = time.time()

  baseTestDir = os.getcwd()

  buildDirName = serialOrMpi+"_"+buildType

  if not performAnyActions(inOptions):
    print "\nNo other actions to perform!\n"
    return success

  print "\nCreating a new build directory if it does not already exist ..."
  createDir(buildDirName)

  absBuildDir = os.path.join(baseTestDir, buildDirName)

  echoChDir(absBuildDir)

  try:

    print ""
    print "A) Get the CMake configure options ("+buildDirName+") ..."
    print ""

    # A.1) Set the base options
  
    cmakeBaseOptions = []
    
    if serialOrMpi == "MPI":
      cmakeBaseOptions.append("-DTPL_ENABLE_MPI:BOOL=ON")
  
    cmakeBaseOptions.append("-DTrilinos_ENABLE_TESTS:BOOL=ON")

    if inOptions.ctestTimeOut:
      cmakeBaseOptions.append(("-DDART_TESTING_TIMEOUT:STRING="+str(inOptions.ctestTimeOut)))
  
    cmakeBaseOptions.extend(extraCMakeOptions)

    readAndAppendCMakeOptions(
      os.path.join("..", getCommonConfigFileName()),
      cmakeBaseOptions)

    readAndAppendCMakeOptions(
      os.path.join("..", getBuildSpecificConfigFileName(serialOrMpi, buildType)),
      cmakeBaseOptions)

    print "\ncmakeBaseOptions:", cmakeBaseOptions

    # A.2) Set the package enable options

    cmakePkgOptions = []
    enablePackagesList = []
  
    if inOptions.enablePackages:
      print "\nEnabling only the explicitly specified packages '"+inOptions.enablePackages+"' ...\n"
      enablePackagesList = inOptions.enablePackages.split(',')
    else:
      diffOutFileName = "../"+getModifiedFilesOutputFileName()
      print "\nDetermining the set of packages to enable by examining "+diffOutFileName+" ...\n"
      if os.path.exists(diffOutFileName):
        updateOutputStr = open(diffOutFileName, 'r').read()
        extractPackageEnablesFromChangeStatus(updateOutputStr, inOptions, enablePackagesList)
      else:
        print "\nThe file "+diffOutFileName+" does not exist!\n"

    for pkg in enablePackagesList:
      cmakePkgOptions.append("-DTrilinos_ENABLE_"+pkg+":BOOL=ON")

    cmakePkgOptions.append("-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON")
  
    if inOptions.enableAllPackages == 'on':
      print "\nEnabling all packages on request ..."
      cmakePkgOptions.append("-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON")
  
    if inOptions.enableFwdPackages:
      print "\nEnabling forward packages on request ..."
      cmakePkgOptions.append("-DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON")
    else:
      cmakePkgOptions.append("-DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=OFF")

    if inOptions.disablePackages:
      print "\nDisabling specified packages '"+inOptions.disablePackages+"' ...\n"
      disablePackagesList = inOptions.disablePackages.split(',')
      for pkg in disablePackagesList:
        cmakePkgOptions.append("-DTrilinos_ENABLE_"+pkg+":BOOL=OFF")

    print "\ncmakePkgOptions:", cmakePkgOptions

    # A.3) Set the combined options

    cmakeOptions    = cmakeBaseOptions + cmakePkgOptions
  
    print "\ncmakeOptions =", cmakeOptions
  
    print "\nCreating base configure file do-configure.base ..."
    createConfigureFile(cmakeBaseOptions, "cmake", inOptions.trilinosSrcDir,
      "do-configure.base")
  
    print "\nCreating package-enabled configure file do-configure ..."
    createConfigureFile(cmakePkgOptions, "./do-configure.base", None, "do-configure")
  
    print ""
    print "B) Do the configuration with CMake ("+buildDirName+") ..."
    print ""

    configurePassed = False
  
    if inOptions.doConfigure:
  
      removeIfExists("CMakeCache.txt")

      cmnd = "./do-configure"
      if inOptions.extraCmakeOptions:
        cmnd += " " + inOptions.extraCmakeOptions

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
    print "C) Do the build ("+buildDirName+") ..."
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

      print "\nFAILED: Skipping the build since configure did not pass!\n"
      
    else:

      print "\nSkipping the build on request!\n"
      if os.path.exists(getBuildSuccessFileName()):
        print "\nA current successful build exists!\n"
        buildPassed = True
      else:
        print "\nFAILED: A current successful build does *not* exist!\n"
  
    print ""
    print "D) Run the tests ("+buildDirName+") ..."
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
        print "\nTest passed!\n"
        echoRunSysCmnd("touch "+getTestSuccessFileName())
      else:
        errStr = "FAILED: ctest failed returning "+str(testRtn)+"!"
        print "\n"+errStr+"\n"
        raise Exception(errStr)
  
    elif inOptions.doTest and not buildPassed:

      print "\nFAILED: Skipping running tests since the build failed!\n"
  
    else:
  
      print "\nSkipping the tests on request!\n"

  except Exception, e:

    success = False
    printStackTrace()

  print ""
  print "E) Analyze the overall results and send email notification ("+buildDirName+") ..."
  print ""

  if performAnyActions(inOptions):

    result = analyzeResultsSendEmail(inOptions, buildDirName,
      enablePackagesList, cmakeOptions, startingTime, timings)
    if not result: succcess = False

  else:

    print "No actions performed, nothing to analyze!"

  return success


def cleanBuildTestCaseOutputFiles(runTestCaseBool, inOptions, baseTestDir, \
  serialOrMpi, buildType ) \
  :

  buildDirName = serialOrMpi+"_"+buildType

  if runTestCaseBool and os.path.exists(buildDirName):

    if inOptions.wipeClean:

      print "\nRemoving the existing build directory "+buildDirName+" (--wipe-clean) ..."
      if os.path.exists(buildDirName):
        echoRunSysCmnd("rm -rf "+buildDirName)

    elif doRemoveOutputFiles(inOptions):

      echoChDir(buildDirName)
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


def runTestCaseDriver(runTestCaseBool, inOptions, baseTestDir, serialOrMpi, buildType,
  extraCMakeOptions, timings ) \
  :

  success = True

  print "\n***"
  print "*** Doing build and test of "+serialOrMpi+" "+buildType+" ..."
  print "***\n"
  
  if runTestCaseBool:

    try:
      echoChDir(baseTestDir)
      writeDefaultBuildSpecificConfigFile(serialOrMpi, buildType)
      result = runTestCase(inOptions, serialOrMpi, buildType,
        extraCMakeOptions, timings.deepCopy())
      if not result: success = False
    except Exception, e:
      success = False
      printStackTrace()

  else:

    print "\nSkipping "+serialOrMpi+" "+buildType+" build/test on request!\n"

  return success

def checkBuildTestCaseStatus(runTestCaseBool, serialOrMpi, buildType, inOptions):

  buildName = serialOrMpi+"_"+buildType

  statusMsg = None

  if not runTestCaseBool:
    buildTestCaseActionsPass = True
    buildTestCaseOkayToCommit = True
    statusMsg = \
      "Test case "+buildName+" was not run! => Does not affect commit/push readiness!"
    return (buildTestCaseActionsPass, buildTestCaseOkayToCommit, statusMsg)

  if not os.path.exists(buildName) and not performAnyBuildTestActions(inOptions):
    buildTestCaseActionsPass = True
    buildTestCaseOkayToCommit = False
    statusMsg = "No configure, build, or test for "+buildName+" was requested!"
    return (buildTestCaseActionsPass, buildTestCaseOkayToCommit, statusMsg)

  if not os.path.exists(buildName):
    buildTestCaseActionsPass = False
    buildTestCaseOkayToCommit = False
    statusMsg = "The directory "+buildName+" does not exist!"

  emailsuccessFileName = buildName+"/"+getEmailSuccessFileName()
  if os.path.exists(emailsuccessFileName):
    buildTestCaseActionsPass = True
  else:
    buildTestCaseActionsPass = False

  testSuccessFileName = buildName+"/"+getTestSuccessFileName()
  if os.path.exists(testSuccessFileName):
    buildTestCaseOkayToCommit = True
  else:
    buildTestCaseOkayToCommit = False

  if not statusMsg:
    statusMsg = getTestCaseSummaryLine(buildName)

  return (buildTestCaseActionsPass, buildTestCaseOkayToCommit, statusMsg)


def getUserCommitMessageStr(inOptions):

  absCommitMsgHeaderFile = inOptions.commitMsgHeaderFile
  if not os.path.isabs(absCommitMsgHeaderFile):
    absCommitMsgHeaderFile = os.path.join(inOptions.trilinosSrcDir, absCommitMsgHeaderFile)

  print "\nExtracting commit message subject and header from the file '" \
        +absCommitMsgHeaderFile+"' ...\n"
  
  commitMsgHeaderFileStr = open(absCommitMsgHeaderFile, 'r').read()
  
  commitEmailBodyStr = commitMsgHeaderFileStr
  
  return commitEmailBodyStr


def getAutomatedStatusSummaryHeaderKeyStr():
  return "Build/Test Cases Summary"


def getAutomatedStatusSummaryHeaderStr(inOptions):
  
  commitEmailBodyStr = "\n\n\n\n" \
      "=============================\n" \
      +getAutomatedStatusSummaryHeaderKeyStr()+"\n" \
      "=============================\n" \
      "\n\n" \
      + getCmndOutput("date", True) + "\n\n" \
  
  return commitEmailBodyStr


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
    if origLogStrList[-3] != "":
      raise Exception("Error, there must be at least one black line before the" \
        " build/test summary block!  This should never happen!")
    origLogStrList = origLogStrList[0:-lastNumBlankLines-1]
  else:
    lastNumBlankLines = -1 # Flag we did not find status header

  return ('\n'.join(origLogStrList), lastNumBlankLines)


def getLastCommitMessageStr(inOptions):

  # Get the raw output from the last current commit log
  rawLogOutput = getCmndOutput(
    "eg cat-file -p HEAD",
    workingDir=inOptions.trilinosSrcDir
    )

  return getLastCommitMessageStrFromRawCommitLogStr(rawLogOutput)[0]


def getLocalCommitsSummariesStr(inOptions):

  # Get the raw output from the last current commit log
  rawLocalCommitsStr = getCmndOutput(
    "eg log --oneline origin..",
    True,
    workingDir=inOptions.trilinosSrcDir
    )

  print \
    "\nLocal commits:" \
    "\n--------------" \
    "\n" \
    + rawLocalCommitsStr

  #print "\nrawLocalCommitsStr:\n=====\n"+rawLocalCommitsStr+"\n=====\n"
  if rawLocalCommitsStr == "\n" or rawLocalCommitsStr == "":
    localCommitsExist = False
  else:
    localCommitsExist = True

  if not localCommitsExist:
    print "\nNo local commits exit!\n"

  localCommitsStr = \
    "\n" \
    "Local commits for this build/test group:\n" \
    "----------------------------------------\n"
  if localCommitsExist:
    localCommitsStr += rawLocalCommitsStr
  else:
    localCommitsStr += "No local commits exist!"
  localCommitsStr += "\n"

  return (localCommitsStr, localCommitsExist)


def checkinTest(inOptions):

  print "\n**********************************************"
  print "*** Performing checkin testing of Trilinos ***"
  print "**********************************************"

  scriptsDir = getScriptBaseDir()
  #print "\nscriptsDir =", scriptsDir

  print "\ntrilinosSrcDir =", inOptions.trilinosSrcDir

  baseTestDir = os.getcwd()
  print "\nbaseTestDir =", baseTestDir

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

  assertEgGitVersions(inOptions)

  success = True

  timings = Timings()

  subjectLine = None

  try:


    print "\n***"
    print "*** 1) Clean old output files ..."
    print "***"
  
    if inOptions.doCommit:
      removeIfExists(getInitialCommitEmailBodyFileName())
      removeIfExists(getInitialCommitOutputFileName())

    if inOptions.doPull:
      removeIfExists(getStatusOutputFileName())
      removeIfExists(getInitialPullOutputFileName())
      removeIfExists(getInitialExtraPullOutputFileName())
      removeIfExists(getInitialPullSuccessFileName())

    removeIfExists(getFinalCommitEmailBodyFileName())
    removeIfExists(getFinalCommitOutputFileName())
    removeIfExists(getCommitStatusEmailBodyFileName())

    removeIfExists(getFinalPullOutputFileName())
    removeIfExists(getModifiedFilesOutputFileName())

    cleanBuildTestCaseOutputFiles(inOptions.withMpiDebug, inOptions, baseTestDir,
      "MPI", "DEBUG" )

    cleanBuildTestCaseOutputFiles(inOptions.withSerialRelease, inOptions, baseTestDir,
      "SERIAL", "RELEASE" )


    print "\n***"
    print "*** 2) Commit changes before pulling updates to merge in ..."
    print "***"

    commitPassed = True # To allow logic below

    if inOptions.doCommit:
    
      try:

        print "\nNOTE: We must commit before doing an 'eg pull --rebase' ...\n"
  
        echoChDir(baseTestDir)

        print "\n2.a) Creating the commit message file ...\n"

        commitEmailBodyStr = getUserCommitMessageStr(inOptions)
        writeStrToFile(getInitialCommitEmailBodyFileName(), commitEmailBodyStr)

        print "\n2.b) Performing the initial local commit (staged and unstaged changes) ...\n"

        print \
          "\nNOTE: This is a trial commit done to allow for a safe pull.  If the  build/test\n"\
          "fails and this commit does not get erased automatically for some reason (i.e.\n" \
          "you killed the script before it was able to finish) then please remove this\n" \
          "commit yourself manually with:\n" \
          "\n" \
          "  $ eg reset --soft HEAD^\n" \
          "\n" \
          "You can then run the checkin-test.py script again to try the test, commit and\n" \
          "push again.\n"

        commitRtn = echoRunSysCmnd(
          "eg commit -a -F "+os.path.join(baseTestDir, getInitialCommitEmailBodyFileName()),
          workingDir=inOptions.trilinosSrcDir,
          outFile=os.path.join(baseTestDir, getInitialCommitOutputFileName()),
          throwExcept=False, timeCmnd=True )

        if commitRtn == 0:
          print "\nCommit passed!\n"
          commitPassed = True
        else:
          print "\nFAILED: Commit failed!\n"
          commitPassed = False

      except Exception, e:
        success = False
        commitPassed = False
        printStackTrace()

    else:

      print "\nSkipping initial commit on request ...\n"

    print "\n***"
    print "*** 3) Update the Trilinos sources ..."
    print "***"

    pullPassed = True

    doingAtLeastOnePull = (inOptions.extraPullFrom or inOptions.doPull)

    if not doingAtLeastOnePull:

      print "\nSkipping all updates on request!\n"

    if inOptions.doCommit and not commitPassed:

      print "\nCommit failed, aborting pull!\n"
      pullPassed = False

    if doingAtLeastOnePull and pullPassed:

      #
      print "\n3.a) Check that there are no uncommited files before doing the pull(s) ...\n"
      #

      statusRtn = echoRunSysCmnd(
        "eg status",
        workingDir=inOptions.trilinosSrcDir,
        outFile=os.path.join(baseTestDir, getStatusOutputFileName()),
        throwExcept=False, timeCmnd=True )
      if statusRtn != 1:
        print "\n'eg status' returned "+str(statusRtn)+": There are uncommitted changes, can not do the pull!\n"
        pullPassed = False

    if doingAtLeastOnePull and pullPassed:

      # NOTE: We want to pull first from the global repo and then from the
      # extra repo so the extra repo's revisions will get rebased on top of
      # the others.  This is what you would want and expect for the remote
      # test/push process where multiple pulls may be needed before it works.

      #
      print "\n3.b) Pull updates from the global 'origin' repo ..."
      #
    
      if inOptions.doPull and pullPassed:
        echoChDir(baseTestDir)
        (updateRtn, updateTimings) = \
          executePull(inOptions, baseTestDir, getInitialPullOutputFileName())
        timings.update += updateTimings
        if updateRtn != 0:
          print "\nPull failed!\n"
          pullPassed = False
      else:
        print "\nSkipping initial pull from 'origin'!\n"

      #
      print "\n3.c) Pull updates from the extra repository '"+inOptions.extraPullFrom+"' ..."
      #

      timings.update = 0
      
      if inOptions.extraPullFrom and pullPassed:
        echoChDir(baseTestDir)
        (updateRtn, updateTimings) = \
          executePull(inOptions, baseTestDir, getInitialExtraPullOutputFileName(),
            inOptions.extraPullFrom )
        timings.update += updateTimings
        if updateRtn != 0:
          print "\nPull failed!\n"
          pullPassed = False
      else:
        print "\nSkipping extra pull from '"+inOptions.extraPullFrom+"'!\n"

    #
    print "\nDetermine overall update pass/fail ...\n"
    #

    echoChDir(baseTestDir)

    # Check for prior successful initial pull
    currentSuccessfullPullExists = os.path.exists(getInitialPullSuccessFileName())

    if (inOptions.doPull or inOptions.extraPullFrom):
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

    if pullPassed:

      getCurrentDiffOutput(inOptions, baseTestDir)

    else:

      print "\nSkipping getting list of modified files because pull failed!\n"


    print "\n***"
    print "*** 5) Running the different build/test cases ..."
    print "***"

    # Determine if we will run the build/test cases or not

    # Set runBuildCases flag
    if not performAnyBuildTestActions(inOptions):
      print "\nNot performing any build cases because no --configure, --build or --test" \
        " was specified!\n"
      runBuildCases = False
    elif inOptions.doCommit and not commitPassed:
      print "\nThe commit failed, skipping running the build/test cases!\n"
      runBuildCases = False
    elif doingAtLeastOnePull:
      if pullPassed:
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

    if runBuildCases:

      echoChDir(baseTestDir)
  
      writeDefaultCommonConfigFile()
  
      result = runTestCaseDriver(inOptions.withMpiDebug, inOptions, baseTestDir,
        "MPI", "DEBUG",
        [
          "-DCMAKE_BUILD_TYPE:STRING=RELEASE",
          "-DTrilinos_ENABLE_DEBUG:BOOL=ON",
          "-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON",
          "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
        ],
        timings
        )
      if not result: success = False
  
      result = runTestCaseDriver(inOptions.withSerialRelease, inOptions, baseTestDir,
        "SERIAL", "RELEASE",
        [
          "-DCMAKE_BUILD_TYPE:STRING=RELEASE",
          "-DTrilinos_ENABLE_DEBUG:BOOL=OFF",
          "-DTrilinos_ENABLE_CHECKED_STL:BOOL=OFF",
          "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF"
        ],
        timings
        )
      if not result: success = False


    print "\n***"
    print "*** 6) Determine overall success and commit/push readiness (and backout commit if failed) ..."
    print "***"

    okayToCommit = False
    okayToPush = False
    forcedCommit = False

    if inOptions.doCommitReadinessCheck or inOptions.doCommit:

      echoChDir(baseTestDir)
  
      okayToCommit = True
      subjectLine = None
      commitEmailBodyExtra = ""

      commitEmailBodyExtra += \
        "\nBuild test results:" \
        "\n-------------------\n"
      buildTestCaseList = [
        [inOptions.withMpiDebug, "MPI", "DEBUG"] ,
        [inOptions.withSerialRelease, "SERIAL", "RELEASE"]
        ]
      for i in range(len(buildTestCaseList)):
        buildTestCase = buildTestCaseList[i]
        (buildTestCaseActionsPass, buildTestCaseOkayToCommit, statusMsg) = \
          checkBuildTestCaseStatus(buildTestCase[0], buildTestCase[1], buildTestCase[2], inOptions)
        print "\n"+statusMsg
        commitEmailBodyExtra += str(i)+") "+buildTestCase[1]+"_"+buildTestCase[2]+" => "+statusMsg
        if not buildTestCaseOkayToCommit:
          commitEmailBodyExtra += " => Not ready for final commit/push!"
        commitEmailBodyExtra += "\n"
        #print "buildTestCaseActionsPass =", buildTestCaseActionsPass
        if not buildTestCaseActionsPass:
          success = False
        if not buildTestCaseOkayToCommit:
          okayToCommit = False
        

      # Print message if a commit is okay or not
      if okayToCommit:
        print "\nThe tests ran and all passed!\n\n" \
          "  => A COMMIT IS OKAY TO BE PERFORMED!"
      else:
        print "\nAt least one of the actions (update, configure, built, test)" \
          " failed or was not performed correctly!\n\n" \
          "  => A COMMIT IS *NOT* OKAY TO BE PERFORMED!"

      # Back out commit if one was performed and buid/test failed
      if not okayToCommit and inOptions.doCommit and commitPassed and not inOptions.forceCommit:
        print "\nNOTICE: Backing out the commit that was just done ...\n"
        try:
          echoRunSysCmnd("eg reset --soft HEAD^",
            workingDir=inOptions.trilinosSrcDir,
            timeCmnd=True,
            throwExcept=False )
        except Exception, e:
          success = False
          printStackTrace()
      
      # Determine if we should do a forced commit/push
      if inOptions.doCommitReadinessCheck and not okayToCommit and commitPassed \
        and inOptions.forceCommit \
        :
        forcedCommitMsg = \
          "\n***" \
          "\n*** WARNING: The acceptance criteria for doing a commit/push has *not*" \
          "\n*** been met, but a commit/push is being forced anyway by --force-commit!" \
          "\n***\n"
        print forcedCommitMsg
        okayToCommit = True
        forcedCommit = True

      # Determine if a push is ready to try or not

      if okayToCommit:
        if currentSuccessfullPullExists:
          print "\nA current successful pull also exists => Ready for final push!\n"
          okayToPush = True
        else:
          commitEmailBodyExtra += \
             "\nA current successful pull does *not* exist => Not ready for final push!\n" \
             "\nExplanation: In order to safely push, the local working directory needs" \
             " to be up-to-date with the global repo or a full integration has not been" \
             " performed!\n"
          print commitEmailBodyExtra
          okayToPush = False
      else:
        okayToPush = False
  
      if okayToPush:
        print "\n  => A PUSH IS READY TO BE PERFORMED!"
      else:
        print "\n  => A PUSH IS *NOT* READY TO BE PERFORMED!"

    else:

      print "\nSkipping commit readiness check on request!"
      okayToCommit = False
      okayToPush = False

  
    print "\n***"
    print "*** 7) Do final push  ..."
    print "***"

    # Attempt the final pull, commit ammend, and push

    pullFinalPassed = True
    ammendFinalCommitPassed = True
    pushPassed = True
    didPush = False
    pushPassed = True
    localCommitSummariesStr = ""

    if not inOptions.doPush:
  
      print "\nNot doing the push on request (--no-push) but sending an email" \
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

        print "\nExplanation: In order to push, the local repo needs to be up-to-date" \
          " with the global repo or the push will not be allowed.  Therefore, a pull" \
          " before the push must be performed if there are updates in the global reop" \
          " regardless if --pull was specified or not.\n"

        (update2Rtn, update2Time) = \
          executePull(inOptions, baseTestDir, getFinalPullOutputFileName())

        if update2Rtn == 0:
          print "\nFinal update passed!\n"
          pullFinalPassed = True
        else:
          print "\nFinal update failed!\n"
          pullFinalPassed = False

        if not pullFinalPassed: okayToPush = False

      #
      print "\n7.b) Ammending the final commit message by appending test results ...\n"
      #

      (localCommitSummariesStr, localCommitsExist) = getLocalCommitsSummariesStr(inOptions)
      #print "\nlocalCommitsExist =", localCommitsExist, "\n"

      if not inOptions.appendTestResults:

        print "\nSkipping appending test results on request (--no-append-test-results)!\n"

      elif not okayToPush:

        print "\nSkippng appending test results due to prior errors!\n"
        ammendFinalCommitPassed = False

      else:  # inOptions.appendTestResults and okayToPush
  
        print "\nAttempting to ammend the final commmit message ...\n"

        try:

          # Get then final commit message
          finalCommitEmailBodyStr = getLastCommitMessageStr(inOptions)
          finalCommitEmailBodyStr += getAutomatedStatusSummaryHeaderStr(inOptions)
          finalCommitEmailBodyStr += commitEmailBodyExtra
          finalCommitEmailBodyStr += localCommitSummariesStr
          if forcedCommit:
            finalCommitEmailBodyStr += (forcedCommitMsg + "\n\n")
          finalCommitEmailBodyStr += getSummaryEmailSectionStr(inOptions)
          writeStrToFile(getFinalCommitEmailBodyFileName(), finalCommitEmailBodyStr)

          # Ammend the final commit message
          if localCommitsExist:
            echoRunSysCmnd(
              "eg commit --amend" \
              " -F "+os.path.join(baseTestDir, getFinalCommitEmailBodyFileName()),
              workingDir=inOptions.trilinosSrcDir,
              outFile=os.path.join(baseTestDir, getFinalCommitOutputFileName()),
              timeCmnd=True
              )
          else:
            print "\nSkipping ammending last commit because there are no local commits!\n"
            
        except Exception, e:
          success = False
          ammendFinalCommitPassed = False
          printStackTrace()

      if not ammendFinalCommitPassed: okayToPush = False

      #
      print "\n7.c) Pushing the the local commits to the global repo ...\n"
      #

      if not okayToPush:

        print "\nNot performing push due to prior errors\n"
        pushPassed = False

      else:
  
        print "\nAttempting to do the push ..."

        pushRtn = echoRunSysCmnd(
          "eg push",
          workingDir=inOptions.trilinosSrcDir,
          outFile=os.path.join(baseTestDir, getPushOutputFileName()),
          throwExcept=False, timeCmnd=True )

        if pushRtn == 0:
          print "\nPush passed!\n"
          pushPassed = True
          didPush = True
        else:
          print "\nPush failed!\n"
          pushPassed = False

      if not pushPassed: okayToPush = False

  
    print "\n***"
    print "*** 8) Create and send push (or readiness status) notification email  ..."
    print "***\n"

    if inOptions.doCommitReadinessCheck:

      #
      print "\n8.a) Getting final status to send out in the summary email ...\n"
      #

      if not commitPassed:
        subjectLine = "INITIAL COMMIT FAILED"
        commitEmailBodyExtra += "\n\nFailed because initial commit failed!" \
          " See '"+getInitialCommitOutputFileName()+"'\n\n"
        success = False
      elif not pullPassed:
        subjectLine = "INITIAL PULL FAILED"
        commitEmailBodyExtra += "\n\nFailed because initial pull failed!" \
          " See '"+getInitialPullOutputFileName()+"'\n\n"
        success = False
      elif not pullFinalPassed:
        subjectLine = "FINAL PULL FAILED"
        commitEmailBodyExtra += "\n\nFailed because the final pull failed!" \
          " See '"+getFinalPullOutputFileName()+"'\n\n"
        success = False
      elif not ammendFinalCommitPassed:
        subjectLine = "AMMEND COMMIT FAILED"
        commitEmailBodyExtra += "\n\nFailed because the final test commit ammend failed!" \
          " See '"+getFinalCommitOutputFileName()+"'\n\n"
        success = False
      elif not pushPassed:
        subjectLine = "PUSH FAILED"
        commitEmailBodyExtra += "\n\nFailed because push failed!" \
          " See '"+getPushOutputFileName()+"'\n\n"
        success = False
      elif inOptions.doPush and pushPassed and forcedCommit:
        subjectLine = "DID FORCED PUSH"
        commitEmailBodyExtra += forcedCommitMsg
        success = True
      elif inOptions.doCommit and commitPassed and forcedCommit:
        subjectLine = "DID FORCED COMMIT"
        commitEmailBodyExtra += forcedCommitMsg
      elif inOptions.doPush:
        if didPush and not forcedCommit:
          subjectLine = "DID PUSH"
        else:
          subjectLine = "ABORTED COMMIT/PUSH"
          commitEmailBodyExtra += "\n\nCommit/push was never attempted since commit/push" \
          " criteria failed!\n\n"
          success = False
      else:
        if okayToPush:
          subjectLine = "READY TO PUSH"
        else:
          subjectLine = "NOT READY TO PUSH"

      #
      print "\n8.b) Create and send out push (or readinessstatus) notification email ..."
      #
  
      subjectLine += ": Trilinos: "+getHostname()
  
      emailBodyStr = subjectLine + "\n\n"
      emailBodyStr += getCmndOutput("date", True) + "\n\n"
      emailBodyStr += commitEmailBodyExtra
      emailBodyStr += localCommitSummariesStr
      emailBodyStr += getSummaryEmailSectionStr(inOptions)
  
      print "\nCommit status email being sent:\n" \
        "--------------------------------\n\n\n\n"+emailBodyStr+"\n\n\n\n"

      summaryCommitEmailBodyFileName = getCommitStatusEmailBodyFileName()
  
      writeStrToFile(summaryCommitEmailBodyFileName, emailBodyStr)
  
      if inOptions.sendEmailTo:
  
        emailAddresses = getEmailAddressesSpaceString(inOptions.sendEmailTo)
        echoRunSysCmnd("sleep 2s")
        echoRunSysCmnd("mailx -s \""+subjectLine+"\" " \
          +emailAddresses+" < "+summaryCommitEmailBodyFileName)
        # Above, we use 'sleep 2s' to try to make sure this email is posted
        # after the last pass/fail email!
  
      else:
  
        print "\nNot sending commit status email because --send-email-to is empty!"

    else:

      print "\nNot performing commit/push or sending out commit/push readiness status on request!"
  
    if not performAnyActions(inOptions) and not inOptions.doPush:

      print \
        "\n***\n" \
        "*** WARNING: No actions were performed!\n" \
        "***\n" \
        "*** Hint: Specify --do-all --commit to perform full integration update/build/test\n" \
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
