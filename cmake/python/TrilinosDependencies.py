#!/bin/env python


import xml.dom.minidom


from GeneralScriptSupport import *


#
# Default file locations
#

defaultTrilinosDepsXmlInFile = getScriptBaseDir()+"/data/TrilinosPackageDependencies.xml"

defaultTrilinosDepsHtmlOutFile = getScriptBaseDir()+"/data/TrilinosPackageDependenciesTable.html"

defaultCDashDepsXmlFile = getScriptBaseDir()+"/data/CDashSubprojectDependencies.xml"


#
# Store and manipulate the dependencies
#


class PackageEmailAddresses:

  def __init__(self, checkin_in, regression_in):
    self.checkin = checkin_in
    self.regression = regression_in

  def __str__(self):
    return "{checkin="+self.checkin+", regression="+self.regression+"}"
  

class PackageDependencies:

  def __init__(self, packageName_in, packageDir_in,
    libRequiredDepPackages_in, libOptionalDepPackages_in,
    testRequiredDepPackages_in, testOptionalDepPackages_in,
    emailAddresses_in, parentPackage_in
    ):
    self.packageName = packageName_in
    self.packageDir = packageDir_in
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

  #print "\n    updatePackageDep("+dep1+", "+dep2+") ..."

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

  #print "\n      newDep =", newDep

  return newDep


class DepStats:
  isDirect = None
  isRequired = None
  isTestDepChain = None
  def __init__(self, isDirect, isRequired, isTestDepChain):
    self.isDirect = isDirect
    self.isRequired = isRequired
    self.isTestDepChain = isTestDepChain


class TrilinosDependencies:


  def __init__(self):
    self.__packagesList = []
    self.__packagesNameToID = {}
    self.__packagesDirToID = {}


  def addPackageDependencies(self, packageDeps):
    packageName = packageDeps.packageName
    packageDir = packageDeps.packageDir
    self.__packagesList.append(packageDeps)
    packageDeps.packageID = len(self.__packagesList)-1 
    self.__packagesNameToID.update( { packageName : packageDeps.packageID } )
    self.__packagesDirToID.update( { packageDir : packageDeps.packageID } )


  def numPackages(self):
    return len(self.__packagesList)


  def packageNameToID(self, packageName):
    return self.__packagesNameToID.get(packageName, -1)


  def getPackageByID(self, packageID):
    return self.__packagesList[packageID]


  def getPackageByName(self, packageName):
    return self.getPackageByID(self.__packagesNameToID[packageName])


  def getPackageByDir(self, packageDir):
    packageID = self.__packagesDirToID.get(packageDir, -1)
    #print "\ngetPackageByDir: packageDir="+packageDir+", packageID="+str(packageID)
    if packageID >= 0:
      return self.__packagesList[packageID]
    return None


  def getPackageNameFromPath(self, fullPath, prefixPath):
    #print "\nfullPath="+fullPath
    fullPathArray = getFilePathArray(fullPath)
    if fullPathArray[0] == "packages":
      regexPathPrefix = "packages/"
      pathPrefix = ""
    else:
      regexPathPrefix = ""
      pathPrefix = "../"
    #print "regexPathPrefix = '"+regexPathPrefix+"'"
    #print "pathPrefix = '"+pathPrefix+"'"
    for packageDep in self.__packagesList:
      regexFilePath = regexPathPrefix+packageDep.packageDir+"/"
      ammendedFullPath = pathPrefix+fullPath 
      #print "\nregexFilePath="+regexFilePath
      #print "ammendedFullPath="+ammendedFullPath
      if re.match(regexFilePath, ammendedFullPath):
        #print "MATCH!"
        return packageDep.packageName
    return u""
    # NOTE: The above loop with match subpackages before it matches
    # packages because subpackages are listed before packages!


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
      #print "\n    updateDepCell: depStats.isDirect="+str(depStats.isDirect)+", depStats.isRequired="+str(depStats.isRequired)+", depCategoryName="+depCategoryName
      newDepName = updatePackageDep(currentDepName, newDepName)

    packageRow[packageID+1] = newDepName


  def updatePackageDepsCategory(self, libsOnly, packageRowID, packageID, depCategory,
    depCategoryName, depStats, trilinosDepsTable
    ):

    packageRow = trilinosDepsTable[packageRowID+1]
    #print "\npackageRow =", packageRow

    depList = getattr(self.__packagesList[packageID], depCategory)
    #print "\ndepList =", depList

    for dep in depList:

      depPackage = self.getPackageByName(dep)
      #print "\n  depPackageName =", depPackage.packageName

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
         trilinosDepsTable)


  def updatePackageDeps(self, libsOnly, packageRowID, packageID, depStats,
    trilinosDepsTable
    ):

    self.updatePackageDepsCategory(libsOnly, packageRowID, packageID,
      "libRequiredDepPackages", "LR", depStats, trilinosDepsTable)
    self.updatePackageDepsCategory(libsOnly, packageRowID, packageID,
      "libOptionalDepPackages", "LO", depStats, trilinosDepsTable)

    # Only process the test dependencies if we are asked to do so
    # (i.e. libsOnly=True) or if this is the top-level package.  The tests for
    # dependent packages are not any kind of dependency for tests for the
    # top-level package.  However, we need to record that these are test
    # dependencies so that any package libraries that get recursed are
    # recorded as 'ITR' or 'ITO' and not as library dependencies.
    if not libsOnly and depStats.isDirect:
      libDepStats = DepStats(True, depStats.isRequired, True)
      self.updatePackageDepsCategory(False, packageRowID, packageID,
        "testRequiredDepPackages", "TR", libDepStats, trilinosDepsTable)
      self.updatePackageDepsCategory(False, packageRowID, packageID,
        "testOptionalDepPackages", "TO", libDepStats, trilinosDepsTable)

  
  def createRawTable(self, libsOnly):

    numPackages = self.numPackages()
    #print "\nnumPackages =", numPackages

    trilinosDepsTable = []

    topRow = [ "Packages" ]
    topRow.extend(["P%02d"%(i+1) for i in range(numPackages)] )
    trilinosDepsTable.append(topRow)

    for packageDeps in self.__packagesList:
      i = packageDeps.packageID
      row = ["P%02d"%(i+1)+") "+packageDeps.packageName]
      row.extend(["" for i in range(numPackages)])
      trilinosDepsTable.append(row)

    for packageDeps in self.__packagesList:
      #print "\npackageName =", packageDeps.packageName
      i = packageDeps.packageID
      trilinosDepsTable[i+1][i+1] = "X"
      self.updatePackageDeps(libsOnly, i, i, DepStats(True, True, False), trilinosDepsTable)

    return trilinosDepsTable

  def createTrilinosPackagesNumberedList(self):
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

    packagesListHtml = self.createTrilinosPackagesNumberedList()

    htmlText = \
      "<p><huge><b>Trilinos Test/Example and Library Package Dependencies</b></huge></p>\n"+\
      "\n"+\
      self.createHtmlFromTable(self.createRawTable(False))+\
      "\n"+\
      packagesListHtml+"\n"+\
      "\n"+\
      "<p><b>Legend</b></p>\n"+\
      "\n"+\
      self.createHtmlTableLegend(False)+\
      "\n"+\
      "<p><b><huge>Trilinos Libary-Only Package Dependencies</huge></b></p>\n"+\
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
      "<title>Trilinos Package Dependencies</title>\n"+\
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


  def writeFullHtmlPage(self, htmlFileName=defaultTrilinosDepsHtmlOutFile):
    htmlString = self.createFullHtmlPage()
    htmlFile = open(htmlFileName, 'w')
    htmlFile.write(htmlString)
    htmlFile.close()


  #
  # CDash stuff
  #


  def createCDashDepsXMLFromRawDepsTable(self, rawTable):
    
    xmlText = ""

    xmlText += "<Project name=\"Trilinos\">\n"

    numPackages = self.numPackages()

    for package_i in range(numPackages):

      packageDeps = self.__packagesList[package_i]

      packageName = packageDeps.packageName

      if packageDeps.parentPackage == "":
        
        xmlText += ("  <SubProject name=\""+packageName+"\">\n")
  
        row = rawTable[package_i+1]
  
        for dep_j in range(numPackages):
          entry = row[dep_j+1]
          if entry and entry != "X":
            depPackageName = self.__packagesList[dep_j].packageName
            depPackageStruct = self.getPackageByName(depPackageName)
            if depPackageStruct.parentPackage == "":
              xmlText += ("    <Dependency name=\""+depPackageName+"\"" + \
                " type=\""+entry+"\"/>\n" )
            else:
              # Don't write subpackage depencencies because CDash should not
              # know about this.  Such dependencies will just not appear but
              # that only affects what CDash displays (and I don't ever look
              # at the subproject dependencies list in CDash really and don't
              # find it very useful given the large number of builds).
              None
  
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

    # rabartl: 2011/07/06: ToDo: Change tha above logic to only
    # write'Dependency' elements for actual Trilinos packages in the place of
    # subpackages.  This will be a little hard to implement and test but we
    # need to do so at some point if we want CDash to know the correct
    # dependencies (but I don't really care).

    return xmlText
  

  def createCDashDepsXML(self):
    return self.createCDashDepsXMLFromRawDepsTable(self.createRawTable(False))


  def writeCDashXmlDepsFile(self, xmlDepsFile=defaultCDashDepsXmlFile):
    xmlString = self.createCDashDepsXML()
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
  checkinEmail = getSingleEmailAddress(emailEle, "Checkin")
  regressionEmail = getSingleEmailAddress(emailEle, "Regression")
  return PackageEmailAddresses(checkinEmail, regressionEmail)


def getParentPackage(packageEle):
  parentPackageEle = packageEle.getElementsByTagName("ParentPackage")[0]
  parentPackage = parentPackageEle.getAttribute('value');
  return parentPackage


def getTrilinosDependenciesFromXmlFile(xmlFile=defaultTrilinosDepsXmlInFile):
  #print "xmlFile =", xmlFile
  packageDepXmlDom = xml.dom.minidom.parse(xmlFile)
  trilinosDependencies = TrilinosDependencies()
  for ele in packageDepXmlDom.childNodes[0].childNodes:
    if ele.nodeType == ele.ELEMENT_NODE:
      packageName = ele.getAttribute('name')
      packageDir = ele.getAttribute('dir')
      #print "\npackageName =", packageName
      packageDeps = PackageDependencies(packageName, packageDir,
        getDependenciesByType(ele, "LIB_REQUIRED_DEP_PACKAGES"),
        getDependenciesByType(ele, "LIB_OPTIONAL_DEP_PACKAGES"),
        getDependenciesByType(ele, "TEST_REQUIRED_DEP_PACKAGES"),
        getDependenciesByType(ele, "TEST_OPTIONAL_DEP_PACKAGES"),
        getPackageEmailAddresses(ele),
        getParentPackage(ele)
        )
      #print "\npackageDeps =", str(packageDeps)
      trilinosDependencies.addPackageDependencies(packageDeps)
  return trilinosDependencies
