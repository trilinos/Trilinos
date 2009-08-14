
#
# General scripting support
#
# NOTE: Included first to check the version of python!
#

#
# ToDo:
#
#  (*) Enable figuring out packages to enable from examining update.out
#  (*) Enable auto-commit


from GeneralScriptSupport import *
import time


def getCommonConfigFileName():
  return "COMMON.config"


def getBuildSpecificConfigFileName(serialOrMpi, buildType):
  return serialOrMpi + "_" + buildType + ".config"


def getUpdateOutputFileName():
  return "update.out"


def getUpdateSuccessFileName():
  return "update.success"


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


def getEmailAddressesSpaceString(emailAddressesCommasStr):
  emailAddressesList = emailAddressesCommasStr.split(',')
  return ' '.join(emailAddressesList)


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
      print "  Appnding: "+line
      cmakeOptions_inout.append(line.strip())


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


def analyzeResultsSendEmail(inOptions, trilinosSrcDir, buildDirName,
  enabledPackagesList, cmakeOptions, startingTime ) \
  :

  print ""
  print "1) Determine what passed and failed ..."
  print ""

  success = False

  # Determine if the configured passed

  configurePassed = None
  configureOutputExists = False

  if inOptions.doConfigure:

    if os.path.exists(getConfigureOutputFileName()):
      configureOutputExists = True

    if os.path.exists(getConfigureSuccessFileName()):
      print "\nThe configure passed!\n"
      configurePassed = True
    else:
      print "\nThe configure FAILED!\n"
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
    else:
      print "\nThe configure FAILED!\n"
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

    else:

      print "\nAt least one of the tests ran FAILED!\n"
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
      subjectLine += ": build FAILED"
      overallPassed = False
      selectedFinalStatus = True

  if inOptions.doConfigure and not selectedFinalStatus:
    if configurePassed:
      subjectLine += ": configure-only passed"
      overallPassed = True
      selectedFinalStatus = True
    elif configureOutputExists:
      subjectLine += ": configure FAILED"
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

  emailBody += "Enabled Packages: " + ', '.join(enabledPackagesList) + "\n"
  emailBody += "Hostname: " + getCmndOutput("hostname", True) + "\n"
  emailBody += "Trilinos Source Dir: " + trilinosSrcDir + "\n"
  emailBody += "Build Dir: " + os.getcwd() + "\n"
  emailBody += "Do Update: " + str(inOptions.doUpdate) + "\n"
  emailBody += "Do Configure: " + str(inOptions.doConfigure) + "\n"
  emailBody += "Do Build: " + str(inOptions.doBuild) + "\n"
  emailBody += "Do Test: " + str(inOptions.doTest) + "\n"
  emailBody += "\nCMake Cache Varibles: " + ' '.join(cmakeOptions) + "\n"
  emailBody += "\n"

  if inOptions.doTest and testOutputExists and numTotalTests:

    emailBody += "Test summary:\n-------------\n\n"

    if inOptions.showAllTests:
      emailBody += getCmndOutput("cat "+getTestOutputFileName())
    else:
      emailBody += getCmndOutput("grep -A 10000 '\% tests passed, ' "+getTestOutputFileName())

  endingTime = time.time()
  totalTime = (endingTime - startingTime) / 60.0

  emailBody += "\n\nFinal:\n------\n\n"
  emailBody += "total time for "+buildDirName+" = "+str(totalTime) + " minutes"

  print "emailBody:\n\n", emailBody

  open(getEmailBodyFileName(),'w').write(emailBody)

  print ""
  print "3) Send the email message ..."
  print ""

  if inOptions.sendEmailTo:

    emailAddresses = getEmailAddressesSpaceString(inOptions.sendEmailTo)
    echoRunSysCmnd("mailx -s \""+subjectLine+"\" "+emailAddresses+" < "+getEmailBodyFileName())

  else:

    print "Not sending email because no email addresses where given!"

  # ToDo: Implement!

  return success


def runTestCase(inOptions, serialOrMpi, buildType, trilinosSrcDir, extraCMakeOptions):

  print "\n***"
  print "*** Doing build and test of "+serialOrMpi+" "+buildType+" ..."
  print "***\n"

  startingTime = time.time()

  baseTestDir = os.getcwd()

  buildDirName = serialOrMpi+"_"+buildType

  if not inOptions.rebuild:
    print "\nRemoving the existing build directory ..."
    if os.path.exists(buildDirName):
      echoRunSysCmnd("rm -rf "+buildDirName) 

  print "Creating a new build directory if it does not already exist ..."
  createDir(buildDirName)

  absBuildDir = os.path.join(baseTestDir, buildDirName)

  echoChDir(absBuildDir)

  removeIfExists(getConfigureSuccessFileName())
  removeIfExists(getBuildSuccessFileName())
  removeIfExists(getTestSuccessFileName())

  try:

    print ""
    print "A) Get the CMake configure options ..."
    print ""

    # A.1) Set the base options
  
    cmakeBaseOptions = [
      "-DCMAKE_BUILD_TYPE:STRING="+buildType
      ]
  
    if serialOrMpi:
  
      cmakeBaseOptions.append("-DTPL_ENABLE_MPI:BOOL=ON")
  
    cmakeBaseOptions.append("-DTrilinos_ENABLE_TESTS:BOOL=ON")
  
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
  
    if inOptions.enablePackages:
      print "\nEnabling the specified packages '"+inOptions.enablePackages+"' ...\n"
      enablePackagesList = inOptions.enablePackages.split(',')
      for pkg in enablePackagesList:
        cmakePkgOptions.append("-DTrilinos_ENABLE_"+pkg+":BOOL=ON")
    else:
      print "\nDetermining the set of packages to enable by examining update.out ...\n"
      # ToDo: Implement this!

    cmakePkgOptions.append("-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON")
  
    if inOptions.enableAllPackages:
      print "\nEnabling all packages on request ..."
      cmakePkgOptions.append("-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON")
  
    if inOptions.enableFwdPackages:
      print "\nEnabling forward packages on request ..."
      cmakePkgOptions.append("-DTrilinos_ENABLE_ALL_FORWARD_DEP_PACAKGES:BOOL=ON")

    print "\ncmakePkgOptions:", cmakePkgOptions

    # A.3) Set the combined options

    cmakeOptions    = cmakeBaseOptions + cmakePkgOptions
  
    print "\ncmakeOptions =", cmakeOptions
  
    print ""
    print "B) Do the configuration with CMake ..."
    print ""
  
    print "Creating base configure file do-configure.base ..."
    createConfigureFile(cmakeBaseOptions, "cmake", trilinosSrcDir, "do-configure.base")
  
    print "Creating package-enabled configure file do-configure ..."
    createConfigureFile(cmakePkgOptions, "./do-configure.base", None, "do-configure")
  
    if inOptions.doConfigure:
  
      removeIfExists("CMakeCache.txt")
  
      echoRunSysCmnd("./do-configure",
        outFile=getConfigureOutputFileName(),
        timeCmnd=True
        )
    
      echoRunSysCmnd("touch "+getConfigureSuccessFileName())
  
    else:
  
      print "\nSkipping configure on request!\n"
  
    print ""
    print "C) Do the build ..."
    print ""
  
    if inOptions.doBuild:
  
      cmnd = "make"
      if inOptions.makeOptions:
        cmnd += " " + inOptions.makeOptions
  
      echoRunSysCmnd(cmnd,
        outFile=getBuildOutputFileName(),
        timeCmnd=True
        )
  
      echoRunSysCmnd("touch "+getBuildSuccessFileName())
  
    else:
  
      print "\nSkipping the build on request ...\n"
  
    print ""
    print "D) Run the tests ..."
    print ""
  
    if inOptions.doTest:
  
      cmnd = "ctest"
      if inOptions.ctestOptions:
        cmnd += " " + inOptions.ctestOptions
  
      echoRunSysCmnd(cmnd,
        outFile=getTestOutputFileName(),
        timeCmnd=True
        )
  
      echoRunSysCmnd("touch "+getTestSuccessFileName())
  
    else:
  
      print "\nSkipping the testing on request!\n"

  except Exception, e:

    traceback.print_exc()

  print ""
  print "E) Analyze the overall results and send email notification ..."
  print ""

  success = analyzeResultsSendEmail(inOptions, trilinosSrcDir, buildDirName,
    enablePackagesList, cmakeOptions, startingTime)

  return success


def runTestCaseDriver(runTestCaseBool, inOptions, baseTestDir, serialOrMpi, buildType,
  trilinosSrcDir, extraCMakeOptions ) \
  :

  success = True
  
  if runTestCaseBool:

    try:
      writeDefaultBuildSpecificConfigFile(serialOrMpi, buildType)
      echoChDir(baseTestDir)
      result = runTestCase(inOptions, serialOrMpi, buildType, trilinosSrcDir, extraCMakeOptions)
      if not result: success = False
    except Exception, e:
      traceback.print_exc()

  else:

    print "\nSkipping "+serialOrMpi+" "+buildType+" build on request!\n"

  return success


def checkinTest(inOptions):

  print "\n***"
  print "*** Performing pre-checkin testing of Trilinos! ..."
  print "***"

  scriptsDir = getScriptBaseDir()
  #print "\nscriptsDir =", scriptsDir

  trilinosSrcDir = '/'.join(scriptsDir.split("/")[0:-2])
  print "\ntrilinosSrcDir =", trilinosSrcDir

  baseTestDir = os.getcwd()
  print "\nbaseTestDir =", baseTestDir

  success = True

  try:

    print "\n***"
    print "*** Update the Trilinos sources ..."
    print "***"
  
    if inOptions.doUpdate:
  
      echoChDir(baseTestDir)
      echoRunSysCmnd(inOptions.updateCommand,
        workingDir=trilinosSrcDir,
        outFile=os.path.join(baseTestDir, "update.out"),
        timeCmnd=True
        )
  
    else:
  
      print "\nSkipping update on request!\n"

    writeDefaultCommonConfigFile()

    result = runTestCaseDriver(inOptions.withMpiDebug, inOptions, baseTestDir,
      "MPI", "DEBUG", trilinosSrcDir,
      [
        "-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON",
        "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
      ]
      )
    if not result: success = False

    result = runTestCaseDriver(inOptions.withSerialRelease, inOptions, baseTestDir,
      "SERIAL", "RELEASE", trilinosSrcDir,
      [
        "-DTrilinos_ENABLE_CHECKED_STL:BOOL=OFF",
        "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF"
      ]
      )
    if not result: success = False

    if inOptions.doCommit:

      if inOptions.doUpdate and inOptions.doConfigure and inOptions.doBuild \
        and inOptions.doTest \
        :

        print "\nToDo: Perform the commit if everything passed!\n"

      else:

        print "\nRefusing to do a commit since --do-update --do-configure" \
          " --do-build and --do-test are not specified!"

        success = false

  except Exception, e:
    traceback.print_exc()

  return success
