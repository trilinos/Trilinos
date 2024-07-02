#!/usr/bin/env python

# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


import xml.dom.minidom
import os
import sys

from FindGeneralScriptSupport import *
from GeneralScriptSupport import *

#
# Default file locations
#


def getDefaultDepsXmlDirectory(rootDir):
  return os.path.join(rootDir, 'cmake', 'dependencies')


def getDefaultDepsXmlInFile(rootDir, projectName):
  return os.path.join(
    getDefaultDepsXmlDirectory(rootDir),
    '%sPackageDependencies.xml' % projectName)


def getDefaultDepsHtmlOutFile(rootDir, projectName):
  return os.path.join(
    getDefaultDepsXmlDirectory(rootDir),
    '%sPackageDependenciesTable.html' % projectName)


def getDefaultCDashDepsXmlFile(rootDir):
  return os.path.join(
    getDefaultDepsXmlDirectory(rootDir),
    'CDashSubprojectDependencies.xml')


#
# Store and manipulate the dependencies
#


class PackageEmailAddresses:

  def __init__(self, regression_in):
    self.regression = regression_in

  def __str__(self):
    return "{regression="+self.regression+"}"
  

class PackageDependencies:

  def __init__(self, packageName_in, packageDir_in, packageType_in,
    libRequiredDepPackages_in, libOptionalDepPackages_in,
    testRequiredDepPackages_in, testOptionalDepPackages_in,
    emailAddresses_in, parentPackage_in
    ):
    self.packageName = packageName_in
    self.packageDir = packageDir_in
    self.packageType = packageType_in
    self.packageID = -1
    self.libRequiredDepPackages = libRequiredDepPackages_in
    self.libOptionalDepPackages = libOptionalDepPackages_in
    self.testRequiredDepPackages = testRequiredDepPackages_in
    self.testOptionalDepPackages = testOptionalDepPackages_in
    self.emailAddresses = emailAddresses_in
    self.parentPackage = parentPackage_in

  def __str__(self):
    return "{\n"+\
      "  packageName="+self.packageName+",\n"+\
      "  packageID="+str(self.packageID)+",\n"+\
      "  packageDir="+str(self.packageDir)+",\n"+\
      "  libRequiredDepPackages="+str(self.libRequiredDepPackages)+",\n"+\
      "  libOptionalDepPackages="+str(self.libOptionalDepPackages)+",\n"+\
      "  testRequiredDepPackages="+str(self.testRequiredDepPackages)+",\n"+\
      "  testOptionalDepPackages="+str(self.testOptionalDepPackages)+" \n"+\
      "  emailAddresses="+str(self.emailAddresses)+"\n"+\
      "  parentPackage="+str(self.parentPackage)+"\n"+\
      "}\n"


def isRequiredDep(dep):
  return (dep[-1] == 'R')


def isDirectDep(dep):
  return (dep[0] != 'I')


def isLibDep(dep):
  return (dep[0] == 'L' or dep[1] == 'L')


#
# (dep1, dep2) => newDep
#
# (*) Required dependencies trump optional dependencies
# (*) Direct dependencies trump indirect dependencies
# (*) Library dependicnes trump test dependencies
#
def updatePackageDep(dep1, dep2):

  #print("\n    updatePackageDep("+dep1+", "+dep2+") ...")

  dep1_required = isRequiredDep(dep1)
  dep1_direct = isDirectDep(dep1)
  dep1_lib = isLibDep(dep1)

  dep2_required = isRequiredDep(dep2)
  dep2_direct = isDirectDep(dep2)
  dep2_lib = isLibDep(dep2)

  selectedDep = False

  if dep1 == dep2:
    newDep = dep1
    selectedDep = True

  # Required trumps optional
  if not selectedDep:
    if dep1_required and not dep2_required:
      newDep = dep1
      selectedDep = True
    elif not dep1_required and dep2_required:
      newDep = dep2
      selectedDep = True

  # Direct trumps indirect
  if not selectedDep:
    if dep1_direct and not dep2_direct:
      newDep = dep1
      selectedDep = True
    elif not dep1_direct and dep2_direct:
      newDep = dep2
      selectedDep = True

  # Library trumps test
  if not selectedDep:
    if dep1_lib and not dep2_lib:
      newDep = dep1
      selectedDep = True
    elif not dep1_lib and dep2_lib:
      newDep = dep2
      selectedDep = True

  assert(selectedDep)

  #print("\n      newDep =", newDep)

  return newDep


class DepStats:
  isDirect = None
  isRequired = None
  isTestDepChain = None
  def __init__(self, isDirect, isRequired, isTestDepChain):
    self.isDirect = isDirect
    self.isRequired = isRequired
    self.isTestDepChain = isTestDepChain


class TribitsDependencies:


  def __init__(self):
    self.__projectName = None
    self.__projectBaseDirName = None
    self.__packagesList = []
    self.__packagesNameToID = {}
    self.__packagesDirToID = {}


  def setProjectName(self, projectName):
    self.__projectName = projectName


  def getProjectName(self):
    return self.__projectName


  def setProjectBaseDirName(self, projectBaseDirName):
    self.__projectBaseDirName = projectBaseDirName


  def getProjectBaseDirName(self):
    return self.__projectBaseDirName


  def addPackageDependencies(self, packageDeps):
    packageName = packageDeps.packageName
    packageDir = packageDeps.packageDir
    self.__packagesList.append(packageDeps)
    packageDeps.packageID = len(self.__packagesList)-1 
    self.__packagesNameToID.update( { packageName : packageDeps.packageID } )
    self.__packagesDirToID.update( { packageDir : packageDeps.packageID } )


  def numPackages(self):
    return len(self.__packagesList)


  def getPackagesNamesList(self, onlyTopLevelPackages=True):
    packagesNamesList = []
    for packageDep in self.__packagesList:
      #print ("packageDep.packageName = "+packageDep.packageName)
      #print ("packageDep.parentPackage = "+packageDep.parentPackage)
      if packageDep.parentPackage == "":
        addPackage = True
      elif not onlyTopLevelPackages:
        addPackage = True
      else:
        addPackage = False
      if addPackage:
        packagesNamesList.append(packageDep.packageName)
    return packagesNamesList


  def packageNameToID(self, packageName):
    return self.__packagesNameToID.get(packageName, -1)


  def getPackageByID(self, packageID):
    return self.__packagesList[packageID]


  def getPackageByName(self, packageName):
    return self.getPackageByID(self.__packagesNameToID[packageName])


  def getPackageByDir(self, packageDir):
    packageID = self.__packagesDirToID.get(packageDir, -1)
    #print("\ngetPackageByDir: packageDir="+packageDir+", packageID="+str(packageID))
    if packageID >= 0:
      return self.__packagesList[packageID]
    return None


  # Note: Path must contain ending "/"
  def getPackageNameFromPath(self, fullPath):
    for packageDep in self.__packagesList:
      regexFilePath = packageDep.packageDir+"/"
      #print("\nregexFilePath="+regexFilePath)
      #print("fullPath="+fullPath)
      if re.match(regexFilePath, fullPath):
        #print("MATCH!")
        return packageDep.packageName
    return u""
    # NOTE: The above loop with match subpackages before it matches
    # packages because subpackages are listed before packages!


  # Returns the paraent package name given a test name
  def getPackageNameFromTestName(self, testName):
    for packageDep in self.__packagesList:
      startTestName = packageDep.packageName+"_"
      #print("\nstartTestName="+startTestName)
      testNameStartIdx = testName.find(startTestName, 0)
      if testNameStartIdx == 0:
        #print("MATCH!")
        if packageDep.parentPackage:
          #print("Subpackage match!")
          return self.getPackageByName(packageDep.parentPackage).packageName
        # Else, is not a subpackage
        return packageDep.packageName
    return u""


  def filterPackageNameList(self, inputPackagesList, keepTypesList, verbose=False):
    if len(inputPackagesList)==1 and inputPackagesList[0]=='':
      return []
    i = 0
    outputPackagesList = []
    for packageName in inputPackagesList:
      #print("packageName = " + packageName)
      if packageName == "ALL_PACKAGES":
        outputPackagesList.append(packageName)
        continue
      packageDep = self.getPackageByName(packageName)
      packageType = packageDep.packageType
      #print("packageType = " + packageType)
      if findInSequence(keepTypesList, packageType) >= 0:
        outputPackagesList.append(packageName)
      else:
        if verbose:
          print(packageName + " of type " + packageType +
                " is being excluded because it is not in the valid list of " +
                "package types [" + ','.join(keepTypesList) + "]")
    return outputPackagesList

  def __str__(self):
    strRep = ""
    for packageDep in self.__packagesList:
      strRep += str(packageDep)
    return strRep


  def updateDepCell(self, packageRow, packageID, depStats, depCategoryName):

    currentDepName = packageRow[packageID+1]

    newDepName = depCategoryName

    # If we are in a test dependency chain, we must change library
    # dependencies to test dependencies.
    if depStats.isTestDepChain:
      newDepName = 'T'+newDepName[1:]

    if depStats.isDirect:
      newDepName = newDepName
    else:
      newDepName = "I"+newDepName

    if not depStats.isRequired:
      newDepName = newDepName[0:-1]+"O"

    if currentDepName:
      #print("\n    updateDepCell: depStats.isDirect="+str(depStats.isDirect)+", depStats.isRequired="+str(depStats.isRequired)+", depCategoryName="+depCategoryName)
      newDepName = updatePackageDep(currentDepName, newDepName)

    packageRow[packageID+1] = newDepName


  def updatePackageDepsCategory(self, libsOnly, packageRowID, packageID, depCategory,
    depCategoryName, depStats, projectDepsTable
    ):

    packageRow = projectDepsTable[packageRowID+1]
    #print("\npackageRow =", packageRow)

    depList = getattr(self.__packagesList[packageID], depCategory)
    #print("\ndepList =", depList)

    for dep in depList:

      depPackage = self.getPackageByName(dep)
      #print("\n  depPackageName =", depPackage.packageName)

      dep_i = depPackage.packageID

      self.updateDepCell(packageRow, dep_i, depStats, depCategoryName)
      
      if not depStats.isRequired:
        isRequiredDep = False
      elif depCategoryName[-1]=="R":
        isRequiredDep = True
      else:
        isRequiredDep = False

      childDepStats = DepStats(False, isRequiredDep, depStats.isTestDepChain)

      self.updatePackageDeps(libsOnly, packageRowID, dep_i, childDepStats,
         projectDepsTable)


  def updatePackageDeps(self, libsOnly, packageRowID, packageID, depStats,
    projectDepsTable
    ):

    self.updatePackageDepsCategory(libsOnly, packageRowID, packageID,
      "libRequiredDepPackages", "LR", depStats, projectDepsTable)
    self.updatePackageDepsCategory(libsOnly, packageRowID, packageID,
      "libOptionalDepPackages", "LO", depStats, projectDepsTable)

    # Only process the test dependencies if we are asked to do so
    # (i.e. libsOnly=True) or if this is the top-level package.  The tests for
    # dependent packages are not any kind of dependency for tests for the
    # top-level package.  However, we need to record that these are test
    # dependencies so that any package libraries that get recursed are
    # recorded as 'ITR' or 'ITO' and not as library dependencies.
    if not libsOnly and depStats.isDirect:
      libDepStats = DepStats(True, depStats.isRequired, True)
      self.updatePackageDepsCategory(False, packageRowID, packageID,
        "testRequiredDepPackages", "TR", libDepStats, projectDepsTable)
      self.updatePackageDepsCategory(False, packageRowID, packageID,
        "testOptionalDepPackages", "TO", libDepStats, projectDepsTable)

  
  def createRawTable(self, libsOnly):

    numPackages = self.numPackages()
    #print("\nnumPackages =", numPackages)

    projectDepsTable = []

    topRow = [ "Packages" ]
    topRow.extend(["P%02d"%(i+1) for i in range(numPackages)] )
    projectDepsTable.append(topRow)

    for packageDeps in self.__packagesList:
      i = packageDeps.packageID
      row = ["P%02d"%(i+1)+") "+packageDeps.packageName]
      row.extend(["" for i in range(numPackages)])
      projectDepsTable.append(row)

    for packageDeps in self.__packagesList:
      #print("\npackageName =", packageDeps.packageName)
      i = packageDeps.packageID
      projectDepsTable[i+1][i+1] = "X"
      self.updatePackageDeps(libsOnly, i, i, DepStats(True, True, False), projectDepsTable)

    return projectDepsTable

  def createProjectPackagesNumberedList(self):
    numPackages = self.numPackages()
    htmlText = "<p><b>Packages:</b> " + \
      ", ".join( \
        [ "P%02d"%(i+1)+") "+self.__packagesList[i].packageName \
           for i in range(self.numPackages())] \
        ) + \
        "</p>"
    return htmlText

  def createHtmlFromTable(self, rawTable):

    numPackages = self.numPackages()

    htmlText = \
      "<TABLE BORDER=4>\n"+\
      "\n"

    for i in range(numPackages+2):
      htmlText += "<COL ALIGN=LEFT>\n"

    topRow = rawTable[0]
    htmlText += "\n<TR>\n"
    for j in range(numPackages+1):
      htmlText += " <TD><b>"+topRow[j]+"</b></TD>\n"
    htmlText += " <TD><b>Packages</b></TD>\n"
    htmlText += "</TR>\n"
      
    for package_i in range(numPackages):
      row = rawTable[package_i+1]
      htmlText += "\n<TR>\n"
      htmlText += " <TD><b>"+row[0]+"</b></TD>\n"
      for j in range(numPackages):
        entry = row[j+1]
        if not entry: entry = "."
        htmlText += " <TD>"+entry+"</TD>\n"
      htmlText += " <TD><b>"+row[0]+"</b></TD>\n"
      htmlText += "</TR>\n"

    htmlText += "\n<TR>\n"
    htmlText += " <TD><b>Packages</b></TD>\n"
    for j in range(numPackages):
      htmlText += " <TD><b>P%02d"%(j+1)+"</b></TD>\n"
    htmlText += " <TD><b>Packages</b></TD>\n"
    htmlText += "</TR>\n"

    htmlText += "</TABLE>\n"
    
    return htmlText


  def createHtmlTableLegend(self, libsOnly):

    htmlText =\
      "\n"+\
      "<ul>\n"+\
      "<li> <b>X</b>: Diagonal entry for the package itself\n"+\
      "<li> <b>LR</b>: Direct library required dependency\n"+\
      "<li> <b>ILR</b>: Indirect library required dependency\n"+\
      "<li> <b>LO</b>: Direct library optional dependency\n"+\
      "<li> <b>ILO</b>: Indirect library optional dependency\n"

    if not libsOnly:
      htmlText +=\
        "<li> <b>TR</b>: Direct test/example required dependency\n"+\
        "<li> <b>ITR</b>: Indirect test/example required dependency\n"+\
        "<li> <b>TO</b>: Direct test/example optional dependency\n"+\
        "<li> <b>ITO</b>: Indirect test/example optional dependency\n"

    htmlText +=\
      "</ul>\n"+\
      "\n"+\
      "NOTE: When more than one type of dependency is present for any cell"+\
      " the final selection is determined in the following order:\n"+\
      "<ul>\n"+\
      "<li> A required dependency trumps an optional dependency\n"+\
      "<li> A direct dependency trumps an indirect dependency\n"+\
      "<li> A library dependency trumps a test/example dependency\n"+\
      "</ul>\n"

    return htmlText


  def createFullHtmlForTables(self):

    packagesListHtml = self.createProjectPackagesNumberedList()

    htmlText = \
      "<p><huge><b>"+self.getProjectName()+" Test/Example and Library Package Dependencies</b></huge></p>\n"+\
      "\n"+\
      self.createHtmlFromTable(self.createRawTable(False))+\
      "\n"+\
      packagesListHtml+"\n"+\
      "\n"+\
      "<p><b>Legend</b></p>\n"+\
      "\n"+\
      self.createHtmlTableLegend(False)+\
      "\n"+\
      "<p><b><huge>"+self.getProjectName()+" Libary-Only Package Dependencies</huge></b></p>\n"+\
      "\n"+\
      self.createHtmlFromTable(self.createRawTable(True))+\
      "\n"+\
      packagesListHtml+"\n"+\
      "\n"+\
      "<p><b>Legend</b></p>\n"+\
      "\n"+\
      self.createHtmlTableLegend(True)

    return htmlText

  def createFullHtmlPage(self):

    htmlText = \
      "<html>\n"+\
      "<head>\n"+\
      "<title>"+self.getProjectName()+" Package Dependencies</title>\n"+\
      "</head>\n"+\
      "\n"+\
      "<body>\n"+\
      "\n"+\
      self.createFullHtmlForTables()+\
      "\n"+\
      "</body>\n"+\
      "\n"+\
      "</html>\n"

    return htmlText


  def writeFullHtmlPage(self, htmlFileName):
    htmlString = self.createFullHtmlPage()
    htmlFile = open(htmlFileName, 'w')
    htmlFile.write(htmlString)
    htmlFile.close()


  #
  # CDash stuff
  #


  def createCDashXML(self):
    
    xmlText = ""

    xmlText += "<Project name=\""+self.getProjectName()+"\">\n"

    projectBaseDirName = self.getProjectBaseDirName()

    numPackages = self.numPackages()

    for package_i in range(numPackages):

      packageDeps = self.__packagesList[package_i]

      packageName = packageDeps.packageName
      packagePath = packageDeps.packageDir

      if packageDeps.parentPackage == "":
        
        xmlText += ("  <SubProject name=\""+packageName+"\">\n")
        
        xmlText += ("    <Path>"+packagePath+"</Path>\n")
  
        xmlText += \
          "    <EmailAddresses>\n"+\
          "      <Email address=\""+packageDeps.emailAddresses.regression+"\"/>\n"+\
          "    </EmailAddresses>\n"

        xmlText += ("  </SubProject>\n")

      else:

        # We don't want to bother CDash with subpackages.  We want CDash testing
        # to operate on the package level!
        None

      # end if

    # end for

    xmlText += "</Project>\n"

    return xmlText


  def writeCDashXmlDepsFile(self, xmlDepsFile):
    xmlString = self.createCDashXML()
    xmlFile = open(xmlDepsFile, 'w')
    xmlFile.write(xmlString)
    xmlFile.close()


#
# Read in the dependencies from XML
#


def getDependenciesByType(packageEle, typeName):
  packageDepsStr = packageEle.getElementsByTagName(typeName)[0].getAttribute('value');
  if len(packageDepsStr) == 0:
    return []
  return packageDepsStr.split(',')


def getSingleEmailAddress(emailEle, emailType):
  singleEmailEle = emailEle.getElementsByTagName(emailType)[0]
  singleEmailAddress = singleEmailEle.getAttribute('address');
  return singleEmailAddress


def getPackageEmailAddresses(packageEle):
  emailEle = packageEle.getElementsByTagName("EmailAddresses")[0]
  regressionEmail = getSingleEmailAddress(emailEle, "Regression")
  return PackageEmailAddresses(regressionEmail)


def getParentPackage(packageEle):
  parentPackageEle = packageEle.getElementsByTagName("ParentPackage")[0]
  parentPackage = parentPackageEle.getAttribute('value');
  return parentPackage


def getProjectDependenciesFromXmlFile(xmlFile):
  packageDepXmlDom = xml.dom.minidom.parse(xmlFile)
  #print("\npackageDepXmlDom =", dir(packageDepXmlDom))
  #print("\npackageDepXmlDom.documentElement =", dir(packageDepXmlDom.documentElement))
  projectDependencies = TribitsDependencies()
  projectDependencies.setProjectName(
    packageDepXmlDom.documentElement.getAttribute('project'))
  projectDependencies.setProjectBaseDirName(
    packageDepXmlDom.documentElement.getAttribute('baseDirName'))
  for ele in packageDepXmlDom.childNodes[0].childNodes:
    if ele.nodeType == ele.ELEMENT_NODE:
      packageName = ele.getAttribute('name')
      packageDir = ele.getAttribute('dir')
      packageType = ele.getAttribute('type')
      #print("\npackageName =", packageName)
      packageDeps = PackageDependencies(packageName, packageDir, packageType,
        getDependenciesByType(ele, "LIB_REQUIRED_DEP_PACKAGES"),
        getDependenciesByType(ele, "LIB_OPTIONAL_DEP_PACKAGES"),
        getDependenciesByType(ele, "TEST_REQUIRED_DEP_PACKAGES"),
        getDependenciesByType(ele, "TEST_OPTIONAL_DEP_PACKAGES"),
        getPackageEmailAddresses(ele),
        getParentPackage(ele)
        )
      #print("\npackageDeps =", str(packageDeps))
      projectDependencies.addPackageDependencies(packageDeps)
  return projectDependencies
