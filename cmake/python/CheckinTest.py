
#
# ToDo:
#
#  (*) Implement checking and control for enableing TrilinosFramework tests or
#  not depending on if anything under cmake/ or packages/teuchos/test/CTestxxx
#  is modified or not, except if TrilinosPackages.cmake and TrilinosTPLs.cmake
#  is specified.
#
#  (*) Implement check that -DTPL_ENABLE* is not specified in any of the
#  standard configure files.  Also, make sure that there are not
#  -DTrilinos_ENABLE commands either.  Force users to do package
#  enables/disables using --enable-packages, --disable-packages.
#
#  (*) Change logic to not enable everything if TrilinosPackages.cmake or
#  TrilinosTPLs.cmake are changed.
#
#  (*) Make this work automatically on branches too.  It should pull and push
#  to 'origin' always but the current branch should always be used.
#
#  (*) Put in checks for the names of Trilinos packages from --enable-packages
#  and --disable-packages arguments.  Right now a mispelled package name would
#  just be ignored.  Also, put in unit tests for this.
#
#  (*) Turn off framework tests by default and turn them in checkin
#      testing ...
#
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
from TrilinosPackageFilePathUtils import *
from GeneralScriptSupport import *
import time


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
  if performAnyBuildTestActions(inOptions) or inOptions.doCommit or inOptions.doPull:
    return True
  return False


def doGenerateOutputFiles(inOptions):
  return performAnyActions(inOptions)


def doRemoveOutputFiles(inOptions):
  return performAnyActions(inOptions)


def executePull(inOptions, baseTestDir, outFile, pullFromRepo=None):
  cmnd = "eg pull --rebase"
  if pullFromRepo:
    print "\nPulling in updates from '"+pullFromRepo+"' ...\n"
    cmnd += " " + pullFromRepo + " master"
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

    writeStrToFile(commonConfigFileStr, commonConfigFileName)


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
      "# NOTE: Please do not add any options here that would select what pacakges\n" \
      "# get enabled or disabled.\n"

    writeStrToFile(buildSpecificConfigFileStr, buildSpecificConfigFileName)


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
    "eg diff --name-status origin/master",
    workingDir=inOptions.trilinosSrcDir,
    outFile=os.path.join(baseTestDir, getModifiedFilesOutputFileName()),
    timeCmnd=True
    )


def extractPackageEnablesFromChangeStatus(updateOutputStr, inOptions_inout,
  enablePackagesList_inout ) \
  :

  trilinosDependencies = getTrilinosDependenciesFromXmlFile(defaultTrilinosDepsXmlInFile)

  modifiedFilesList = extractFilesListMatchingPattern(
    updateOutputStr.split('\n'), reModifedFiles )

  for modifiedFileFullPath in modifiedFilesList:

    if isGlobalBuildFile(modifiedFileFullPath):
      if inOptions_inout.enableAllPackages == 'default':
        print "\nModifed file: '"+modifiedFileFullPath+"'\n" \
          "  => Enabling all Trilinos packages!"
        inOptions_inout.enableAllPackages = 'on'

    packageName = getPackageNameFromPath(trilinosDependencies, modifiedFileFullPath)
    if packageName and findInSequence(enablePackagesList_inout, packageName) == -1:
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
  
    writeStrToFile(doConfigStr, configFileName)
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
      print "\nThe tests when never even run!\n"
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
    emailBody += "\nExtra CMake Options: " + inOptions.extraCmakeOptions + "\n"
  if inOptions.makeOptions:
    emailBody += "\nMake Options: " + inOptions.makeOptions + "\n"
  if inOptions.ctestOptions:
    emailBody += "\nCTest Options: " + inOptions.ctestOptions + "\n"
  emailBody += "\n"
  emailBody += getStageStatus("Update", inOptions.doPull, updatePassed, timings.update)
  emailBody += getStageStatus("Configure", inOptions.doConfigure, configurePassed, timings.configure)
  emailBody += getStageStatus("Build", inOptions.doBuild, buildPassed, timings.build)
  emailBody += getStageStatus("Test", inOptions.doTest, testsPassed, timings.test)
  emailBody += "\n"

  if inOptions.doTest and testOutputExists and numTotalTests:

    emailBody += "Test summary:\n-------------\n\n"

    if inOptions.showAllTests:
      emailBody += getCmndOutput("cat "+getTestOutputFileName())
    else:
      emailBody += getCmndOutput("grep -A 10000 '\% tests passed, ' "+getTestOutputFileName())

  else:

    emailBody += "\n***\n*** WARNING: There are no test results!\n***\n\n"

  endingTime = time.time()
  totalTime = (endingTime - startingTime) / 60.0

  emailBody += "\n\nFinal:\n------\n\n"
  emailBody += "Total time for "+buildDirName+" = "+str(totalTime) + " minutes"

  #print "emailBody:\n\n\n\n", emailBody, "\n\n\n\n"

  writeStrToFile(emailBody, getEmailBodyFileName())

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


def getTestCaseEmailSummary(doTestCaseBool, testCaseName):
  summaryEmailSectionStr = ""
  if doTestCaseBool:
    summaryEmailSectionStr += \
      "\n\n"+testCaseName+" Results:\n" \
      "------------------------\n" \
      "\n"
    absEmailBodyFileName = testCaseName+"/"+getEmailBodyFileName()
    if os.path.exists(absEmailBodyFileName):
      testCaseEmailStrArray = open(absEmailBodyFileName, 'r').readlines()
      for line in testCaseEmailStrArray:
        summaryEmailSectionStr += "  " + line
      summaryEmailSectionStr += "\n"
    else:
        summaryEmailSectionStr += \
          "Error, The build/test was never completed!" \
          " (the file '"+absEmailBodyFileName+"' does not exist.)\n"
  return summaryEmailSectionStr


def getSummaryEmailSectionStr(inOptions):
  summaryEmailSectionStr = \
    "\n\n-------------------\n" \
    "Summary of Results:\n" \
    "-------------------\n"
  summaryEmailSectionStr += \
    getTestCaseEmailSummary(inOptions.withMpiDebug, "MPI_DEBUG")
  summaryEmailSectionStr += \
    getTestCaseEmailSummary(inOptions.withSerialRelease, "SERIAL_RELEASE")
  return summaryEmailSectionStr

  
def runTestCase(inOptions, serialOrMpi, buildType, extraCMakeOptions,
  timings ) \
  :

  success = True

  startingTime = time.time()

  baseTestDir = os.getcwd()

  buildDirName = serialOrMpi+"_"+buildType

  if not inOptions.rebuild:
    print "\nRemoving the existing build directory ..."
    if os.path.exists(buildDirName):
      echoRunSysCmnd("rm -rf "+buildDirName)

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
  
      print "\nSkipping the testing on request!\n"

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


def cleanTestCaseOutputFiles(runTestCaseBool, inOptions, baseTestDir, \
  serialOrMpi, buildType ) \
  :
  buildDirName = serialOrMpi+"_"+buildType
  if runTestCaseBool and doRemoveOutputFiles(inOptions) \
    and os.path.exists(buildDirName) \
    :
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


def checkBuildCheckinStatus(runTestCaseBool, serialOrMpi, buildType):

  buildName = serialOrMpi+"_"+buildType

  statusMsg = ""

  if not runTestCaseBool:
    return (True,
      "\nTest case "+buildName+" was not run!  Does not affect commit/push readiness!\n")

  if not os.path.exists(buildName):
    statusMsg += "\nThe directory "+buildName+" does not exist!  Not ready for final commit/push!\n"
    return (False, statusMsg)

  testSuccessFileName = buildName+"/"+getTestSuccessFileName()
  if not os.path.exists(testSuccessFileName):
     statusMsg += "\nThe file "+testSuccessFileName+" does not exist!  Not ready for final commit/push!\n"
     return (False, statusMsg)

  emailSuccessFileName = buildName+"/"+getEmailSuccessFileName()
  if not os.path.exists(emailSuccessFileName):
    statusMsg += "\nThe file "+emailSuccessFileName+" does not exist!  Not ready for final commit!\n"
    return (False, statusMsg)

  statusMsg += "\nThe tests successfully passed for "+buildName+"!\n"
  return (True, statusMsg)


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
  return "Automated status information"


def getAutomatedStatusSummaryHeaderStr(inOptions):
  
  commitEmailBodyStr = "\n\n\n\n" \
      "=============================\n" \
      +getAutomatedStatusSummaryHeaderKeyStr()+"\n" \
      "=============================\n" \
      "\n\n" \
      + getCmndOutput("date", True) + "\n\n" \
  
  return commitEmailBodyStr


def getLastCommitMessageStr(inOptions):

  # Get the raw output from the last current commit log
  rawLogOutput = getCmndOutput(
    "eg cat-file -p HEAD",
    workingDir=inOptions.trilinosSrcDir
    )

  # Extract the original log message
  origLogStrList = []
  pastHeader = False
  numBlankLines = 0
  foundStatusHeader = False
  for line in rawLogOutput.split('\n'):
    #print "\nline = '"+line+"'\n"
    if pastHeader:
      origLogStrList.append(line)
      if line == "":
        numBlankLines += 1
      else:
        numBlankLines = 0
      if line == getAutomatedStatusSummaryHeaderKeyStr():
        foundStatusHeader = True
        break
    if line == "":
      pastHeader = True

  if foundStatusHeader:
    origLogStrList = origLogStrList[0:-numBlankLines-2]

  return '\n'.join(origLogStrList)


def getLocalCommitsSummariesStr(inOptions):

  # Get the raw output from the last current commit log
  rawLocalCommitsStr = getCmndOutput(
    " eg log --oneline origin/master..master",
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
    "\n\n\n\n" \
    "----------------------------------------\n" \
    "Local commits for this build/test group:\n" \
    "----------------------------------------\n\n"
  if localCommitsExist:
    localCommitsStr += rawLocalCommitsStr
  else:
    localCommitsStr += "No local commits exist!"

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
      removeIfExists(getInitialPullSuccessFileName())

    removeIfExists(getFinalCommitEmailBodyFileName())
    removeIfExists(getFinalCommitOutputFileName())
    removeIfExists(getCommitStatusEmailBodyFileName())

    removeIfExists(getFinalPullOutputFileName())
    removeIfExists(getModifiedFilesOutputFileName())

    cleanTestCaseOutputFiles(inOptions.withMpiDebug, inOptions, baseTestDir,
      "MPI", "DEBUG" )

    cleanTestCaseOutputFiles(inOptions.withSerialRelease, inOptions, baseTestDir,
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
        writeStrToFile(commitEmailBodyStr, getInitialCommitEmailBodyFileName())

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
      print "\n3.d) Determine overall update pass/success ...\n"
      #

      echoChDir(baseTestDir)

      if (inOptions.doPull or inOptions.extraPullFrom):
        if pullPassed:
          print "\nUpdate passed!\n"
          echoRunSysCmnd("touch "+getInitialPullSuccessFileName())
        else:
          print "\nUpdate failed!\n"
      elif os.path.exists(getInitialPullSuccessFileName()):
        print "\nA previous update was performed and was successful!"
        pullPassed = True
      elif inOptions.allowNoPull:
        print "\nNot performing update since --skip-update was passed in\n"
        pullPassed = True
      else:
        print "\nNo previous successful update is still current!"
        pullPassed = False


    print "\n***"
    print "*** 4) Get the list of all the modified files ..."
    print "***"

    getCurrentDiffOutput(inOptions, baseTestDir)


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
    print "*** 6) Determine overall commit/push readiness ..."
    print "***"

    if inOptions.doCommitReadinessCheck or inOptions.doCommit:

      echoChDir(baseTestDir)
  
      okToCommit = True
      subjectLine = None
      commitEmailBodyExtra = ""
  
      (buildOkay, statusMsg) = \
        checkBuildCheckinStatus(inOptions.withMpiDebug, "MPI", "DEBUG")
      print statusMsg
      commitEmailBodyExtra += statusMsg
      if not buildOkay:
        okToCommit = False
        
      (buildOkay, statusMsg) = \
        checkBuildCheckinStatus(inOptions.withSerialRelease, "SERIAL", "RELEASE")
      print statusMsg
      commitEmailBodyExtra += statusMsg
      if not buildOkay:
        okToCommit = False
  
  
      if okToCommit:
        print "\nThe tests ran and all passed!\n\n" \
          "  => A PUSH IS OKAY TO BE PERFORMED!"
      else:
        print "\nAt least one of the actions (update, configure, built, test)" \
          " failed or was not performed correctly!\n\n" \
          "  => A PUSH IS *NOT* READY TO BE PERFORMED!"

    else:

      print "\nSkipping commit readiness check on request!"
      okToCommit = False
  
    print "\n***"
    print "*** 7) Do commit and push  ..."
    print "***"

    # Back out commit if one was performed and buid/test failed
    if not okToCommit and inOptions.doCommit and commitPassed and not inOptions.forceCommit:
      print "\nNOTICE: Backing out the commit that was just done ...\n"
      try:
        echoRunSysCmnd("eg reset --soft HEAD^",
          workingDir=inOptions.trilinosSrcDir,
          timeCmnd=True,
          throwExcept=False )
      except Exception, e:
        success = False
        printStackTrace()
    
    # Determine if we should do a forced commit
    forcedCommit = False
    if inOptions.doCommitReadinessCheck and not okToCommit and commitPassed \
      and inOptions.forceCommit \
      :
      forcedCommitMsg = \
        "\n***" \
        "\n*** WARNING: The acceptance criteria for doing a commit/push has *not*" \
        "\n*** been met, but a commit/push is being forced anyway by --force-commit!" \
        "\n***\n"
      print forcedCommitMsg
      okToCommit = True
      forcedCommit = True

    # Attempt the final pull, commit ammend, and push

    pullFinalPassed = True
    ammendFinalCommitPassed = True
    pushPassed = True
    didPush = False
    localCommitSummariesStr = ""
    okToPush = okToCommit

    if inOptions.doCommitReadinessCheck and okToCommit:

      #
      print "\n7.a) Performing a final pull to make sure there are no conflicts ...\n"
      #
      
      if not inOptions.doPull:

        print "\nSkipping the final pull (--skip-final-pull)!\n"

      elif not okToPush:

        print "\nSkippng final pull due to prior errors!\n"
        pullFinalPassed = False

      else: # inOptions.doPull and okToPush

        (update2Rtn, update2Time) = \
          executePull(inOptions, baseTestDir, getFinalPullOutputFileName())

        if update2Rtn == 0:
          print "\nFinal update passed!\n"
          pullFinalPassed = True
        else:
          print "\nFinal update failed!\n"
          pullFinalPassed = False

      if not pullFinalPassed: okToPush = False

      #
      print "\n7.b) Ammending the final commit message by appending test results ...\n"
      #

      (localCommitSummariesStr, localCommitsExist) = getLocalCommitsSummariesStr(inOptions)
      #print "\nlocalCommitsExist =", localCommitsExist, "\n"

      if not inOptions.appendTestResults:

        print "\nSkipping apending test results on request (--no-append-test-results)!\n"

      elif not okToPush:

        print "\nSkippng apending test results due to prior errors!\n"
        ammendFinalCommitPassed = False

      else:  # inOptions.appendTestResults and okToPush
  
        print "\nAttempting to ammend the final commmit message ...\n"

        try:

          # Get then final commit message
          finalCommitEmailBodyStr = getLastCommitMessageStr(inOptions)
          finalCommitEmailBodyStr += getAutomatedStatusSummaryHeaderStr(inOptions)
          finalCommitEmailBodyStr += localCommitSummariesStr
          if forcedCommit:
            finalCommitEmailBodyStr += (forcedCommitMsg + "\n\n")
          finalCommitEmailBodyStr += getSummaryEmailSectionStr(inOptions)
          writeStrToFile(finalCommitEmailBodyStr, getFinalCommitEmailBodyFileName())

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

      if not ammendFinalCommitPassed: okToPush = False

      #
      print "\n7.c) Pushing the the local commits to the global repo ...\n"
      #

      pushPassed = True

      if not inOptions.doPush:
  
        print "\nNot doing the push on request (--no-push) but sending an email" \
              " about the commit/push readiness status ..."
  
        if okToCommit:
          subjectLine = "READY TO PUSH"
        else:
          subjectLine = "NOT READY TO PUSH"

      elif not okToPush:

        print "\nNot performing push due to prior errors\n"
        pushPassed = False

      else: # inOptions.doPush and okToPush:
  
        print "\nAttempting to do a push ..."

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

      if not pushPassed: okToPush = False

    else:

      print "\nNot attempted final commit and/or push!"

  
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
      elif not pullPassed:
        subjectLine = "INITIAL PULL FAILED"
        commitEmailBodyExtra += "\n\nFailed because initial pull failed!" \
          " See '"+getInitialPullOutputFileName()+"'\n\n"
      elif not pullFinalPassed:
        subjectLine = "FINAL PULL FAILED"
        commitEmailBodyExtra += "\n\nFailed because the final pull failed!" \
          " See '"+getFinalPullOutputFileName()+"'\n\n"
      elif not ammendFinalCommitPassed:
        subjectLine = "AMMEND COMMIT FAILED"
        commitEmailBodyExtra += "\n\nFailed because the final test commit ammend failed!" \
          " See '"+getFinalCommitOutputFileName()+"'\n\n"
      elif not pushPassed:
        subjectLine = "PUSH FAILED"
        commitEmailBodyExtra += "\n\nFailed because push failed!" \
          " See '"+getPushOutputFileName()+"'\n\n"
      elif inOptions.doPush and pushPassed and forcedCommit:
        subjectLine = "FORCED COMMIT/PUSH"
        commitEmailBodyExtra += forcedCommitMsg
      elif inOptions.doCommit and commitPassed and forcedCommit:
        subjectLine = "FORCED COMMIT"
        commitEmailBodyExtra += forcedCommitMsg
      elif inOptions.doPush:
        if didPush and not forcedCommit:
          subjectLine = "DID PUSH"
        else:
          subjectLine = "ABORTED COMMIT/PUSH"
          commitEmailBodyExtra += "\n\nCommit/push was never attempted since commit/push" \
          " criteria failed!\n\n"
      else:
        if okToCommit:
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
  
      writeStrToFile(emailBodyStr, summaryCommitEmailBodyFileName)
  
      if inOptions.sendEmailTo:
  
        emailAddresses = getEmailAddressesSpaceString(inOptions.sendEmailTo)
        echoRunSysCmnd("sleep 2s && mailx -s \""+subjectLine+"\" " \
          +emailAddresses+" < "+summaryCommitEmailBodyFileName)
        # Above, we use 'sleep 2s' to try to make sure this email is posted
        # after the last pass/fail email!
  
      else:
  
        print "\nNot sending commit status email because --send-email-to is empty!"

    else:

      print "\nNot performing commit/push or sending out commit/push readiness status on request!"
  
    if not performAnyActions(inOptions) and not inOptions.doCommit:

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

  # Print the final status at the very end
  if subjectLine:
    print \
      "\n\n" + subjectLine + "\n\n"

  return success



#
# Unit testing code
# 


import unittest



class MockOptions:
  def __init__(self):
    self.enableAllPackages = 'default'


trilinosDependencies = getTrilinosDependenciesFromXmlFile(defaultTrilinosDepsXmlInFile)


class testCheckinTest(unittest.TestCase):


  def test_isGlobalBuildFile_00(self):
    self.assertEqual( isGlobalBuildFile( 'Trilinos_version.h' ), True )


  def test_isGlobalBuildFile_01(self):
    self.assertEqual( isGlobalBuildFile( 'CMakeLists.txt' ), True )


  def test_isGlobalBuildFile_02(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/TrilinosPackages.cmake' ), True )


  def test_isGlobalBuildFile_03(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/TrilinosCMakeQuickstart.txt' ), False )


  def test_isGlobalBuildFile_04(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/ctest/experimental_build_test.cmake' ),
      False )


  def test_isGlobalBuildFile_05(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/DependencyUnitTests/blah' ),
      False )


  def test_isGlobalBuildFile_06(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/TPLs/FindTPLBLAS.cmake' ),
      True )


  def test_isGlobalBuildFile_07(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/TPLs/FindTPLLAPACK.cmake' ),
      True )


  def test_isGlobalBuildFile_08(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/TPLs/FindTPLMPI.cmake' ),
      True )


  def test_isGlobalBuildFile_09(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/TPLs/FindTPLDummy.cmake' ),
      False )


  def test_isGlobalBuildFile_10(self):
    self.assertEqual( isGlobalBuildFile( 'cmake/utils/SetNotFound.cmake' ),
      True )


  def test_extractPackageEnablesFromChangeStatus_1(self):

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
      enablePackagesList)

    self.assertEqual( options.enableAllPackages, 'on' )
    self.assertEqual( enablePackagesList, ['Teuchos'] )


  def test_extractPackageEnablesFromChangeStatus_2(self):

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
      enablePackagesList)

    self.assertEqual( options.enableAllPackages, 'default' )
    self.assertEqual( enablePackagesList, ['NOX', 'Thyra'] )


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(testCheckinTest))
    return suite


if __name__ == '__main__':
  unittest.TextTestRunner(verbosity=2).run(suite())
