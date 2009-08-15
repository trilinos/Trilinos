
#
# ToDo:
#
#  (*) Remove every output file/directory first with --from-scrarch ...
#  (*) Enable auto-commit with embedded email.out files at end ...
#  (*) Make it clear that all code and tests and examples should pass before committing!
#  (*) Turn off framework tests by default and turn them in pre-checkin testing ...
#  (*) Turn off generation of HTML/XML files by default and turn them on in pre-checkin testing ...
#

#
# General scripting support
#
# NOTE: Included first to check the version of python!
#

from TrilinosDependencies import getTrilinosDependenciesFromXmlFile
from TrilinosDependencies import defaultTrilinosDepsXmlInFile
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
      print "  Appnding: "+line.strip()
      cmakeOptions_inout.append(line.strip())


reModifedFiles = re.compile(r"^M (.+)$")


def isGlobalCmakeBuildFile(modifiedFileFullPathArray):
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


def getPackageNameFromPathArray(trilinosDependencies, modifiedFileFullPathArray):
  if modifiedFileFullPathArray[0] == "packages":
    packageDir = modifiedFileFullPathArray[1]
    packageStruct = trilinosDependencies.getPackageByDir(packageDir)
    if packageStruct:
      return packageStruct.packageName
  return ""


def extractPackageEnablesFromChangeStatus(updateOutputStr, inOptions_inout,
   enablePackagesList_inout )\
  :

  trilinosDependencies = getTrilinosDependenciesFromXmlFile(defaultTrilinosDepsXmlInFile)

  for updateLine in updateOutputStr.split('\n'):
    
    #print "\nupdateLine =", updateLine

    reModifiedFilesMatch = reModifedFiles.match(updateLine)

    if reModifiedFilesMatch:

      modifedFileFullPath = reModifiedFilesMatch.group(1).strip()
      #print "\nmodifedFileFullPath =", modifedFileFullPath

      modifiedFileFullPathArray = modifedFileFullPath.split('/')
      #print "\nmodifiedFileFullPathArray =", modifiedFileFullPathArray

      if isGlobalCmakeBuildFile(modifiedFileFullPathArray):

        if inOptions_inout.enableAllPackages == 'default':
          print "\nModifed file: '"+modifedFileFullPath+"'\n" \
            "  => Enabling all Trilinos packages!"
          inOptions_inout.enableAllPackages = 'on'

      packageName = getPackageNameFromPathArray(trilinosDependencies, modifiedFileFullPathArray)
      if packageName and findInSequence(enablePackagesList_inout, packageName) == -1:
        print "\nModified file: '"+modifedFileFullPath+"'\n" \
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
  if inOptions.extraCmakeOptions:
    emailBody += "\nExtra CMake Options: " + inOptions.extraCmakeOptions + "\n"
  if inOptions.ctestOptions:
    emailBody += "\nCTest Options: " + inOptions.ctestOptions + "\n"
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
  emailBody += "Total time for "+buildDirName+" = "+str(totalTime) + " minutes"

  print "emailBody:\n\n\n\n", emailBody, "\n\n\n\n"

  open(getEmailBodyFileName(),'w').write(emailBody)

  print ""
  print "3) Send the email message ..."
  print ""

  if inOptions.sendEmailTo:

    emailAddresses = getEmailAddressesSpaceString(inOptions.sendEmailTo)
    echoRunSysCmnd("mailx -s \""+subjectLine+"\" "+emailAddresses+" < "+getEmailBodyFileName())

  else:

    print "Not sending email because no email addresses where given!"

  # 3) Return final result

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
      "-DCMAKE_BUILD_TYPE:STRING="+buildType,
      "-DTrilinos_ALLOW_NO_PACKAGES:BOOL=OFF"
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
    enablePackagesList = []
  
    if inOptions.enablePackages:
      print "\nEnabling the specified packages '"+inOptions.enablePackages+"' ...\n"
      enablePackagesList = inOptions.enablePackages.split(',')
    else:
      print "\nDetermining the set of packages to enable by examining update.out ...\n"
      updateOutputStr = open("../update.out", 'r').read()
      extractPackageEnablesFromChangeStatus(updateOutputStr, inOptions, enablePackagesList)

    for pkg in enablePackagesList:
      cmakePkgOptions.append("-DTrilinos_ENABLE_"+pkg+":BOOL=ON")

    cmakePkgOptions.append("-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON")
  
    if inOptions.enableAllPackages == 'on':
      print "\nEnabling all packages on request ..."
      cmakePkgOptions.append("-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON")
  
    if inOptions.enableFwdPackages:
      print "\nEnabling forward packages on request ..."
      cmakePkgOptions.append("-DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON")

    if inOptions.extraCmakeOptions:
      print "\nAppending extra CMake options from command-line ..."
      cmakePkgOptions.extend(inOptions.extraCmakeOptions.split(" "))

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


  return (allKeywords, ' '. join(testDirs))


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
    self.assertEqual( isGlobalCmakeBuildFile( ['CMakeLists.txt'] ), True )


  def test_isGlobalCmakeBuildFile_02(self):
    self.assertEqual( isGlobalCmakeBuildFile( ['cmake', 'TrilinosPackages.cmake' ] ), True )


  def test_isGlobalCmakeBuildFile_03(self):
    self.assertEqual( isGlobalCmakeBuildFile( ['cmake', 'TrilinosCMakeQuickstart.txt' ] ), False )


  def test_isGlobalCmakeBuildFile_04(self):
    self.assertEqual( isGlobalCmakeBuildFile( ['cmake', 'ctest', 'experimental_build_test.cmake' ] ),
      False )


  def test_isGlobalCmakeBuildFile_05(self):
    self.assertEqual( isGlobalCmakeBuildFile( ['cmake', 'DependencyUnitTests', 'blah' ] ),
      False )


  def test_isGlobalCmakeBuildFile_06(self):
    self.assertEqual( isGlobalCmakeBuildFile( ['cmake', 'TPLs', 'FindTPLBLAS.cmake' ] ),
      True )


  def test_isGlobalCmakeBuildFile_07(self):
    self.assertEqual( isGlobalCmakeBuildFile( ['cmake', 'TPLs', 'FindTPLLAPACK.cmake' ] ),
      True )


  def test_isGlobalCmakeBuildFile_08(self):
    self.assertEqual( isGlobalCmakeBuildFile( ['cmake', 'TPLs', 'FindTPLMPI.cmake' ] ),
      True )


  def test_isGlobalCmakeBuildFile_09(self):
    self.assertEqual( isGlobalCmakeBuildFile( ['cmake', 'TPLs', 'FindTPLDummy.cmake' ] ),
      False )


  def test_isGlobalCmakeBuildFile_10(self):
    self.assertEqual( isGlobalCmakeBuildFile( ['cmake', 'utils', 'SetNotFound.cmake' ] ),
      True )


  def test_getPackageNameFromPathArray_01(self):
    self.assertEqual(
      getPackageNameFromPathArray( trilinosDependencies, ['packages', 'teuchos', 'CMakeLists.txt' ] ), 'Teuchos' )


  def test_getPackageNameFromPathArray_02(self):
    self.assertEqual(
      getPackageNameFromPathArray( trilinosDependencies, ['packages', 'thyra', 'src', 'blob.cpp' ] ), 'Thyra' )


  def test_getPackageNameFromPathArray_03(self):
    self.assertEqual(
      getPackageNameFromPathArray( trilinosDependencies, ['packages', 'blob', 'blob' ] ), '' )


  def test_extractPackageEnablesFromChangeStatus_1(self):

    updateOutputStr = """
? packages/tpetra/doc/html
? packages/trilinoscouplings/doc/html
? packages/triutils/doc/html
? sampleScripts/checkin-test-gabriel.sh
M cmake/TrilinosPackages.cmake
M cmake/python/checkin-test.py
M cmake/python/dump-cdash-deps-xml-file.py
M packages/teuchos/example/ExplicitInstantiation/four_files/CMakeLists.txt
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
M packages/nox/src/dummy.C
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
