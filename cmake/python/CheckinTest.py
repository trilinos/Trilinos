
#
# ToDo:
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


def getTestCaseNamee(serialOrMpi, buildType):
  return serialOrMpi + "_" + buildType


def getBuildSpecificConfigFileName(serialOrMpi, buildType):
  return getTestCaseNamee(serialOrMpi, buildType) + ".config"


def getUpdateOutputFileName():
  return "update.out"


def getUpdateSuccessFileName():
  return "update.success"


def getUpdateOutput2FileName():
  return "update2.out"


def getConfigureOutputFileName():
  return "do-configure.out"


def getConfigureSuccessFileName():
  return "do-configure.success"


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


def getCommitEmailBodyFileName():
  return "commitEmailBody.out"


def getCommitStatusEmailBodyFileName():
  return "commitStatusEmailBody.out"


def getCommitOutputFileName():
  return "commit.out"


def getHostname():
  return getCmndOutput("hostname", True)


def getEmailAddressesSpaceString(emailAddressesCommasStr):
  emailAddressesList = emailAddressesCommasStr.split(',')
  return ' '.join(emailAddressesList)


def performAnyActions(inOptions):
  if inOptions.doUpdate or inOptions.doConfigure or inOptions.doBuild \
    or inOptions.doTest or inOptions.doAll \
    :
    return True
  return False


def doGenerateOutputFiles(inOptions):
  return performAnyActions(inOptions)


def doRemoveOutputFiles(inOptions):
  return performAnyActions(inOptions)


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

    open(commonConfigFileName, 'w').write(commonConfigFileStr)


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

    open(buildSpecificConfigFileName, 'w').write(buildSpecificConfigFileStr)


def readAndAppendCMakeOptions(fileName, cmakeOptions_inout):

  if not os.path.exists(fileName):
    return

  print "\nAppending options from "+fileName+":"

  cmakeOptionsFile = open(fileName, 'r')

  for line in cmakeOptionsFile:
    if line[0] != '#':
      print "  Appnding: "+line.strip()
      cmakeOptions_inout.append(line.strip())


reModifedFiles = re.compile(r"^[MA] (.+)$")


def isGlobalCmakeBuildFile(modifiedFileFullPath):
  modifiedFileFullPathArray = getFilePathArray(modifiedFileFullPath)
  if len(modifiedFileFullPathArray)==1 and modifiedFileFullPathArray[0] == "CMakeLists.txt":
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


def extractPackageEnablesFromChangeStatus(updateOutputStr, inOptions_inout,
  enablePackagesList_inout ) \
  :

  trilinosDependencies = getTrilinosDependenciesFromXmlFile(defaultTrilinosDepsXmlInFile)

  modifiedFilesList = extractFilesListMatchingPattern(
    updateOutputStr.split('\n'), reModifedFiles )

  for modifiedFileFullPath in modifiedFilesList:

    if isGlobalCmakeBuildFile(modifiedFileFullPath):
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
  
    open(configFileName, 'w').write(doConfigStr)
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


def analyzeResultsSendEmail(inOptions, trilinosSrcDir, buildDirName,
  enabledPackagesList, cmakeOptions, startingTime, timings ) \
  :

  print ""
  print "1) Determine what passed and failed ..."
  print ""

  success = False

  # Determine if the update passed

  updatePassed = None
  updateOutputExists = False

  if inOptions.doUpdate:

    if os.path.exists("../"+getUpdateOutputFileName()):
      updateOutputExists = True

    if os.path.exists("../"+getUpdateSuccessFileName()):
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

  if inOptions.doTest:
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

  if inOptions.doUpdate and not selectedFinalStatus:
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
    subjectLine = "passed : " + subjectLine
  else:
    subjectLine = "FAILED : " + subjectLine

  print "\nsubjectLine = '"+subjectLine+"'\n"

  success = overallPassed

  # 2.b) Construct the email body

  emailBody = subjectLine + "\n\n"

  emailBody += getCmndOutput("date", True) + "\n\n"

  emailBody += "Enabled Packages: " + ', '.join(enabledPackagesList) + "\n"
  if inOptions.disablePackages:
    emailBody += "Disabled Packages: " + inOptions.disablePackages + "\n"
  emailBody += "Hostname: " + getHostname() + "\n"
  emailBody += "Source Dir: " + trilinosSrcDir + "\n"
  emailBody += "Build Dir: " + os.getcwd() + "\n"
  emailBody += "\nCMake Cache Varibles: " + ' '.join(cmakeOptions) + "\n"
  if inOptions.extraCmakeOptions:
    emailBody += "\nExtra CMake Options: " + inOptions.extraCmakeOptions + "\n"
  if inOptions.makeOptions:
    emailBody += "\nMake Options: " + inOptions.makeOptions + "\n"
  if inOptions.ctestOptions:
    emailBody += "\nCTest Options: " + inOptions.ctestOptions + "\n"
  emailBody += "\n"
  emailBody += getStageStatus("Update", inOptions.doUpdate, updatePassed, timings.update)
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

  open(getEmailBodyFileName(),'w').write(emailBody)

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
        summaryEmailSectionStr += "Error, the file '"+absEmailBodyFileName+"' does not exist!\n"
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

  
def runTestCase(inOptions, serialOrMpi, buildType, trilinosSrcDir, extraCMakeOptions,
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
    print "A) Get the CMake configure options ..."
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
      print "\nEnabling the specified packages '"+inOptions.enablePackages+"' ...\n"
      enablePackagesList = inOptions.enablePackages.split(',')
    else:
      print "\nDetermining the set of packages to enable by examining update.out ...\n"
      updateOutFileName = "../"+getUpdateOutputFileName()
      if os.path.exists(updateOutFileName):
        updateOutputStr = open(updateOutFileName, 'r').read()
        extractPackageEnablesFromChangeStatus(updateOutputStr, inOptions, enablePackagesList)
      else:
        print "\nThe file "+updateOutFileName+" does not exist!\n"

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
    createConfigureFile(cmakeBaseOptions, "cmake", trilinosSrcDir, "do-configure.base")
  
    print "\nCreating package-enabled configure file do-configure ..."
    createConfigureFile(cmakePkgOptions, "./do-configure.base", None, "do-configure")
  
    print ""
    print "B) Do the configuration with CMake ..."
    print ""
  
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
      else:
        print "\nConfigure failed returning "+str(configureRtn)+"!\n"
        raise Exception("Configure failed!")
  
    else:
  
      print "\nSkipping configure on request!\n"
  
    print ""
    print "C) Do the build ..."
    print ""
  
    if inOptions.doBuild:
  
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
      else:
        print "\nBuild failed returning "+str(buildRtn)+"!\n"
        raise Exception("Build failed!")
  
    else:

      print "\nSkipping the build on request!\n"
  
    print ""
    print "D) Run the tests ..."
    print ""
  
    if inOptions.doTest:
  
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
        print "\nTest failed returning "+str(testRtn)+"!\n"
        raise Exception("Test failed!")
  
    else:
  
      print "\nSkipping the testing on request!\n"

  except Exception, e:

    success = False

    traceback.print_exc()

  print ""
  print "E) Analyze the overall results and send email notification ..."
  print ""

  if performAnyActions(inOptions):

    result = analyzeResultsSendEmail(inOptions, trilinosSrcDir, buildDirName,
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
    if inOptions.doConfigure or inOptions.doUpdate:
      removeIfExists(getConfigureOutputFileName())
      removeIfExists(getConfigureSuccessFileName())
    if inOptions.doBuild or inOptions.doConfigure or inOptions.doUpdate:
      removeIfExists(getBuildOutputFileName())
      removeIfExists(getBuildSuccessFileName())
    if inOptions.doTest or inOptions.doBuild or inOptions.doConfigure or inOptions.doUpdate:
      removeIfExists(getTestOutputFileName())
      removeIfExists(getTestSuccessFileName())
    removeIfExists(getEmailBodyFileName())
    removeIfExists(getEmailSuccessFileName())
    echoChDir("..")


def runTestCaseDriver(runTestCaseBool, inOptions, baseTestDir, serialOrMpi, buildType,
  trilinosSrcDir, extraCMakeOptions, timings ) \
  :

  success = True

  print "\n***"
  print "*** Doing build and test of "+serialOrMpi+" "+buildType+" ..."
  print "***\n"
  
  if runTestCaseBool:

    try:
      echoChDir(baseTestDir)
      writeDefaultBuildSpecificConfigFile(serialOrMpi, buildType)
      result = runTestCase(inOptions, serialOrMpi, buildType, trilinosSrcDir,
        extraCMakeOptions, timings.deepCopy())
      if not result: success = False
    except Exception, e:
      success = False
      traceback.print_exc()

  else:

    print "\nSkipping "+serialOrMpi+" "+buildType+" build/test on request!\n"

  return success


def checkBuildCheckinStatus(runTestCaseBool, serialOrMpi, buildType):

  buildName = serialOrMpi+"_"+buildType

  statusMsg = ""

  if not runTestCaseBool:
    return (True,
      "\nTest case "+buildName+" was not run!  Does not affect commit readiness!\n")

  if not os.path.exists(buildName):
    statusMsg += "\nThe directory "+buildName+" does not exist!  Not ready to commit!\n"
    return (False, statusMsg)

  testSuccessFileName = buildName+"/"+getTestSuccessFileName()
  if not os.path.exists(testSuccessFileName):
     statusMsg += "\nThe file "+testSuccessFileName+" does not exist!  Not ready to commit!\n"
     return (False, statusMsg)

  emailSuccessFileName = buildName+"/"+getEmailSuccessFileName()
  if not os.path.exists(emailSuccessFileName):
    statusMsg += "\nThe file "+emailSuccessFileName+" does not exist!  Not ready to commit!\n"
    return (False, statusMsg)

  statusMsg += "\nThe tests successfully passed for "+buildName+"!\n"
  return (True, statusMsg)


def checkinTest(inOptions):

  print "\n**********************************************"
  print "*** Performing checkin testing of Trilinos ***"
  print "**********************************************"

  scriptsDir = getScriptBaseDir()
  #print "\nscriptsDir =", scriptsDir

  trilinosSrcDir = '/'.join(scriptsDir.split("/")[0:-2])
  print "\ntrilinosSrcDir =", trilinosSrcDir

  baseTestDir = os.getcwd()
  print "\nbaseTestDir =", baseTestDir

  if inOptions.doAll:
    inOptions.doUpdate = True
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

    removeIfExists(getCommitEmailBodyFileName())
    removeIfExists(getCommitStatusEmailBodyFileName())
    removeIfExists(getCommitOutputFileName())

    if inOptions.doUpdate:
      removeIfExists(getUpdateOutputFileName())
      removeIfExists(getUpdateSuccessFileName())

    removeIfExists(getUpdateOutput2FileName())

    cleanTestCaseOutputFiles(inOptions.withMpiDebug, inOptions, baseTestDir,
      "MPI", "DEBUG" )

    cleanTestCaseOutputFiles(inOptions.withSerialRelease, inOptions, baseTestDir,
      "SERIAL", "RELEASE" )


    print "\n***"
    print "*** 2) Update the Trilinos sources ..."
    print "***"
  
    if inOptions.doUpdate:
    
      try:
  
        echoChDir(baseTestDir)
        (updateRtn, timings.update) = echoRunSysCmnd(inOptions.updateCommand,
          workingDir=trilinosSrcDir,
          outFile=os.path.join(baseTestDir, getUpdateOutputFileName()),
          timeCmnd=True, returnTimeCmnd=True, throwExcept=False
          )

        if updateRtn == 0:
          "\nUpdate passed!\n"
          echoRunSysCmnd("touch "+getUpdateSuccessFileName())
          updatePassed = True
        else:
          "\nUpdate failed!\n"
          updatePassed = False

      except Exception, e:
        success = False
	updatePassed = False
        traceback.print_exc()
  
    else:
  
      print "\nSkipping update on request!\n"

      if os.path.exists(getUpdateSuccessFileName()):
        print "\nA previous update was performed and was successful!"
        updatePassed = True
      else:
        print "\nA previous update was *not* performed or was *not* successful!"
        updatePassed = False


    print "\n***"
    print "*** 3) Running the different build cases ..."
    print "***"

    if updatePassed:

      echoChDir(baseTestDir)
  
      writeDefaultCommonConfigFile()
  
      result = runTestCaseDriver(inOptions.withMpiDebug, inOptions, baseTestDir,
        "MPI", "DEBUG", trilinosSrcDir,
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
        "SERIAL", "RELEASE", trilinosSrcDir,
        [
          "-DCMAKE_BUILD_TYPE:STRING=RELEASE",
          "-DTrilinos_ENABLE_DEBUG:BOOL=OFF",
          "-DTrilinos_ENABLE_CHECKED_STL:BOOL=OFF",
          "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF"
        ],
        timings
        )
      if not result: success = False

    else:

      print "\nNot doing any builds because the update was not performed or failed!\n"

    print "\n***"
    print "*** 4) Determine overall commit readiness ..."
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
          "  => A COMMIT IS OKAY TO BE PERFORMED!"
      else:
        print "\nAt least one of the actions (update, configure, built, test)" \
          " failed or was not performed correctly!\n\n" \
          "  => A COMMIT IS *NOT* READY TO BE PERFORMED!"

    else:

      print "\nSkipping commit readiness check on request!"
      okToCommit = False
  
    print "\n***"
    print "*** 5) Do commit or send email about commit readiness status ..."
    print "***"

    if inOptions.doCommitReadinessCheck:

      if inOptions.doCommit:
  
        print "\nAttempting to do a commit ..."
  
        forcedCommit = False
  
        if not okToCommit and inOptions.forceCommit:
          forcedCommitMsg = \
            "\n***" \
            "\n*** WARNING: The acceptance criteria for doing a commit has *not*" \
            "\n*** been met, but a commit is being forced anyway!" \
            "\n***\n"
          print forcedCommitMsg
          okToCommit = True
          forcedCommit = True
  
        if okToCommit:
  
          print "\nDoing a last update to avoid not-up-to-date status ...\n"
  
          if inOptions.doFinalUpdate:
  
            update2Rtn = echoRunSysCmnd(inOptions.updateCommand,
              workingDir=trilinosSrcDir,
              outFile=os.path.join(baseTestDir, getUpdateOutput2FileName()),
              throwExcept=False,
              timeCmnd=True
              )
            if update2Rtn != 0: okToCommit = False
  
          else:
            
            print "\nSkipping the final update on request!\n"
            update2Rtn = 0
  
          absCommitMsgHeaderFile = inOptions.commitMsgHeaderFile
          if not os.path.isabs(absCommitMsgHeaderFile):
            absCommitMsgHeaderFile = os.path.join(trilinosSrcDir,absCommitMsgHeaderFile)
  
          print "\nExtracting commit message subject and header from the file '" \
            +absCommitMsgHeaderFile+"' ...\n"
  
          commitMsgHeaderFileStr = open(absCommitMsgHeaderFile, 'r').read()
  
          commitEmailBodyStr = commitMsgHeaderFileStr
  
          commitEmailBodyStr += "\n\n\n\n" \
            "=============================\n" \
            "Automated status information\n" \
            "=============================\n" \
            "\n\n" \
            + getCmndOutput("date", True) + "\n\n" \
  
          if forcedCommit:
            commitEmailBodyStr += (forcedCommitMsg + "\n\n")
  
          commitEmailBodyStr += \
            getSummaryEmailSectionStr(inOptions)
  
          commitMsgFile = getCommitEmailBodyFileName()
          open(commitMsgFile, 'w').write(commitEmailBodyStr)
  
          if okToCommit:
            commitRtn = echoRunSysCmnd(
              "cvs commit -F "+os.path.join(baseTestDir,commitMsgFile),
              workingDir=trilinosSrcDir,
              outFile=os.path.join(baseTestDir,getCommitOutputFileName()),
              throwExcept=False,
              timeCmnd=True
              )
  
          if update2Rtn != 0:
            okToCommit = False
            subjectLine = "COMMIT FAILED"
            commitEmailBodyExtra += "\n\nCommit failed because final update failed!  See 'update2.out'\n\n"
          elif commitRtn == 0:
            if forcedCommit:
              subjectLine = "FORCED COMMIT"
              commitEmailBodyExtra += forcedCommitMsg
            else:
              subjectLine = "DID COMMIT"
          else:
            subjectLine = "COMMIT FAILED"
            commitEmailBodyExtra += "\n\nCommit failed!  See the file 'commit.out'\n\n"
  
        else:
  
          subjectLine = "ABORTED COMMIT"
  
          commitEmailBodyExtra += "\n\nCommit was never attempted since commit criteria failed!\n\n"
  
      else:
  
        print "\nNot doing the commit but sending an email about the commit readiness status ..."
  
        if okToCommit:
          subjectLine = "READY TO COMMIT"
        else:
          subjectLine = "NOT READY TO COMMIT"
  
      print "\nCreate and send out commit (readiness) status notification email ..."
  
      subjectLine += ": Trilinos: "+getHostname()
  
      if not updatePassed:
        commitEmailBodyExtra += "The update failed!  See the file 'update.out'!\n"
  
      emailBodyStr = subjectLine + "\n\n"
      emailBodyStr += getCmndOutput("date", True) + "\n\n"
      emailBodyStr += commitEmailBodyExtra
      emailBodyStr += getSummaryEmailSectionStr(inOptions)
  
      print "\nCommit status email being sent:\n" \
        "--------------------------------\n\n\n\n"+emailBodyStr+"\n\n\n\n"
  
      summaryCommitEmailBodyFileName = getCommitStatusEmailBodyFileName()
      open(summaryCommitEmailBodyFileName, 'w').write(emailBodyStr)
  
      if inOptions.sendEmailTo:
  
        emailAddresses = getEmailAddressesSpaceString(inOptions.sendEmailTo)
        echoRunSysCmnd("sleep 2s && mailx -s \""+subjectLine+"\" " \
          +emailAddresses+" < "+summaryCommitEmailBodyFileName)
        # Above, we use 'sleep 2s' to try to make sure this email is posted
        # after the last pass/fail email!
  
      else:
  
        print "\nNot sending commit status email because --send-email-to is empty!"

    else:

      print "\nSkipping commit or sending commit readiness status on request!"
  
    if not performAnyActions(inOptions) and not inOptions.doCommit:

      print "\n***"
      print "*** WARNING: No actions were performed!"
      print "***"
      print "*** Specify --do-all to perform full test or --commit to commit a previously run test!"
      print "***\n"
  
  except Exception, e:

    success = False

    traceback.print_exc()

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


  def test_isGlobalCmakeBuildFile_01(self):
    self.assertEqual( isGlobalCmakeBuildFile( 'CMakeLists.txt' ), True )


  def test_isGlobalCmakeBuildFile_02(self):
    self.assertEqual( isGlobalCmakeBuildFile( 'cmake/TrilinosPackages.cmake' ), True )


  def test_isGlobalCmakeBuildFile_03(self):
    self.assertEqual( isGlobalCmakeBuildFile( 'cmake/TrilinosCMakeQuickstart.txt' ), False )


  def test_isGlobalCmakeBuildFile_04(self):
    self.assertEqual( isGlobalCmakeBuildFile( 'cmake/ctest/experimental_build_test.cmake' ),
      False )


  def test_isGlobalCmakeBuildFile_05(self):
    self.assertEqual( isGlobalCmakeBuildFile( 'cmake/DependencyUnitTests/blah' ),
      False )


  def test_isGlobalCmakeBuildFile_06(self):
    self.assertEqual( isGlobalCmakeBuildFile( 'cmake/TPLs/FindTPLBLAS.cmake' ),
      True )


  def test_isGlobalCmakeBuildFile_07(self):
    self.assertEqual( isGlobalCmakeBuildFile( 'cmake/TPLs/FindTPLLAPACK.cmake' ),
      True )


  def test_isGlobalCmakeBuildFile_08(self):
    self.assertEqual( isGlobalCmakeBuildFile( 'cmake/TPLs/FindTPLMPI.cmake' ),
      True )


  def test_isGlobalCmakeBuildFile_09(self):
    self.assertEqual( isGlobalCmakeBuildFile( 'cmake/TPLs/FindTPLDummy.cmake' ),
      False )


  def test_isGlobalCmakeBuildFile_10(self):
    self.assertEqual( isGlobalCmakeBuildFile( 'cmake/utils/SetNotFound.cmake' ),
      True )


  def test_extractPackageEnablesFromChangeStatus_1(self):

    updateOutputStr = """
? packages/tpetra/doc/html
? packages/trilinoscouplings/doc/html
? packages/triutils/doc/html
? sampleScripts/checkin-test-gabriel.sh
M cmake/TrilinosPackages.cmake
M cmake/python/checkin-test.py
M cmake/python/dump-cdash-deps-xml-file.py
P packages/thyra/dummy.blah
A packages/teuchos/example/ExplicitInstantiation/four_files/CMakeLists.txt
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
M cmake/python/checkin-test.py
M cmake/python/dump-cdash-deps-xml-file.py
A packages/nox/src/dummy.C
P packages/stratimikos/dummy.blah
M packages/thyra/src/Thyra_ConfigDefs.hpp
M packages/thyra/CMakeLists.txt
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
