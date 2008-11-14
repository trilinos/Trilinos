#!/bin/env python


import xml.dom.minidom


from GeneralScriptSupport import *


#
# Default file locations
#

defaultTrilinosDepsXmlInFile = getScriptsDir()+"/data/TrilinosPackageDependencies.xml"

defaultTrilinosDepsHtmlOutFile = getScriptsDir()+"/data/TrilinosPackageDependenciesTable.html"


#
# Store and manipulate the dependencies
#

class PackageDependencies:

  packageName = None
  packageID = None
  libRequiredDepPackages = None
  libOptionalDepPackages = None
  testRequiredDepPackages = None
  testOptionalDepPackages = None

  def __init__(self, packageName_in,
    libRequiredDepPackages_in, libOptionalDepPackages_in,
    testRequiredDepPackages_in, testOptionalDepPackages_in
    ):
    self.packageName = packageName_in
    self.packageID = -1
    self.libRequiredDepPackages = libRequiredDepPackages_in
    self.libOptionalDepPackages = libOptionalDepPackages_in
    self.testRequiredDepPackages = testRequiredDepPackages_in
    self.testOptionalDepPackages = testOptionalDepPackages_in

  def __str__(self):
    return "{\n"+\
      "  packageName="+self.packageName+",\n"+\
      "  packageID="+str(self.packageID)+",\n"+\
      "  libRequiredDepPackages="+str(self.libRequiredDepPackages)+",\n"+\
      "  libOptionalDepPackages="+str(self.libOptionalDepPackages)+",\n"+\
      "  testRequiredDepPackages="+str(self.testRequiredDepPackages)+",\n"+\
      "  testOptionalDepPackages="+str(self.testOptionalDepPackages)+" \n"+\
      "}\n"


class TrilinosDependencies:


  __packagesList = None
  __packagesDirToList = None


  def __init__(self):
    self.__packagesList = []
    self.__packagesDirToList = {}


  def addPackageDependencies(self, packageDeps):
    packageName = packageDeps.packageName
    self.__packagesList.append(packageDeps)
    packageDeps.packageID = len(self.__packagesList)-1 
    self.__packagesDirToList.update(      { packageName : packageDeps.packageID } )


  def numPackages(self):
    return len(self.__packagesList)


  def getPackageByName(self, packageName):
    return self.__packagesList[self.__packagesDirToList[packageName]]


  def __str__(self):
    strRep = ""
    for packageDep in self.__packagesList:
      strRep += str(packageDep)
    return strRep


  def updateDepCell(self, packageRow, packageID, isDirect, isRequired, depCategoryName):

    currentDepName = packageRow[packageID+1]

    if currentDepName and currentDepName[-1] == 'R':

      newDepName = currentDepName

    else:
  
      if isDirect:
        newDepName = depCategoryName
      else:
        newDepName = "I"+depCategoryName
  
      if not isRequired:
        newDepName = newDepName[0:-1]+"O"

    packageRow[packageID+1] = newDepName


  def updatePackageDepsCategory(self, libsOnly, packageRowID, packageID, depCategory,
    depCategoryName, isDirect, isRequired, trilinosDepsTable
    ):
    packageRow = trilinosDepsTable[packageRowID+1]
    #print "\npackageRow =", packageRow
    depList = getattr(self.__packagesList[packageID], depCategory)
    #print "\ndepList =", depList
    for dep in depList:
      depPackage = self.getPackageByName(dep)
      #print "\ndepPackage =", depPackage
      dep_i = depPackage.packageID
      self.updateDepCell(packageRow, dep_i, isDirect, isRequired, depCategoryName)
      #packageRow[dep_i+1] = depCategoryName # ToDo: Use a function to do this!
      if depCategoryName[-1]=="R":
        isRequiredDep = True
      else:
        isRequiredDep = False
      self.updatePackageDeps(libsOnly, packageRowID, dep_i, False, isRequiredDep,
         trilinosDepsTable)


  def updatePackageDeps(self, libsOnly, packageRowID, packageID, isDirect, isRequired,
    trilinosDepsTable
    ):
    self.updatePackageDepsCategory(libsOnly, packageRowID, packageID,
      "libRequiredDepPackages", "LR", isDirect, isRequired, trilinosDepsTable)
    self.updatePackageDepsCategory(libsOnly, packageRowID, packageID,
      "libOptionalDepPackages", "LO", isDirect, isRequired, trilinosDepsTable)
    if not libsOnly:
      self.updatePackageDepsCategory(False, packageRowID, packageID,
        "testRequiredDepPackages", "TR", isDirect, isRequired, trilinosDepsTable)
      self.updatePackageDepsCategory(False, packageRowID, packageID,
        "testOptionalDepPackages", "TO", isDirect, isRequired, trilinosDepsTable)

  
  def createRawTable(self, libsOnly):

    numPackages = self.numPackages()
    #print "\nnumPackages =", numPackages

    trilinosDepsTable = []

    topRow = [ "Packages" ]
    topRow.extend(["P"+str(i+1) for i in range(numPackages)] )
    trilinosDepsTable.append(topRow)

    for packageDeps in self.__packagesList:
      i = packageDeps.packageID
      row = ["P"+str(i+1)+") "+packageDeps.packageName]
      row.extend(["" for i in range(numPackages)])
      trilinosDepsTable.append(row)

    for packageDeps in self.__packagesList:
      i = packageDeps.packageID
      trilinosDepsTable[i+1][i+1] = "X"
      self.updatePackageDeps(libsOnly, i, i, True, True, trilinosDepsTable)

    return trilinosDepsTable


  def createHtmlFromTable(self, rawTable):

    numPackages = self.numPackages()

    htmlText = \
      "<TABLE BORDER="+str(numPackages)+">\n"+\
      "\n"

    for i in range(numPackages+1):
      htmlText += "<COL ALIGN=LEFT>\n"

    topRow = rawTable[0]
    htmlText += "\n<TR>\n"
    for j in range(numPackages+1):
      htmlText += " <TD><b>"+topRow[j]+"</b><TD>\n"
    htmlText += "</TR>\n"
      
    for package_i in range(numPackages):
      row = rawTable[package_i+1]
      htmlText += "\n<TR>\n"
      htmlText += " <TD><b>"+row[0]+"</b><TD>\n"
      for j in range(numPackages):
        entry = row[j+1]
        if not entry: entry = "."
        htmlText += " <TD>"+entry+"<TD>\n"
      htmlText += "</TR>\n"

    htmlText += "</TABLE>\n"
    
    return htmlText


  def createHtmlTableLegend(self, libsOnly):

    htmlText =\
      "\n"+\
      "<ul>\n"+\
      "\n"+\
      "<li> <b>X</b>: Diagonal entry for the package itself\n"+\
      "<li> <b>LR</b>: Direct library required dependency\n"+\
      "<li> <b>ILR</b>: Indirect library required dependency\n"+\
      "<li> <b>LO</b>: Dirrect library optional dependency\n"+\
      "<li> <b>ILO</b>: Indirect library optional dependency\n"

    if not libsOnly:
      htmlText +=\
        "<li> <b>TR</b>: Direct test/example required dependency\n"+\
        "<li> <b>ITR</b>: Indirect test/example required dependency\n"+\
        "<li> <b>TO</b>: Dirrect test/example optional dependency\n"+\
        "<li> <b>ITO</b>: Indirect test/example optional dependency\n"

    htmlText +=\
      "\n"+\
      "</ul>\n"

    return htmlText


  def createFullHtmlForTables(self):

    htmlText = \
      "<p><b>Trilinos Libary Package Dependencies</b></p>\n"+\
      "\n"+\
      self.createHtmlFromTable(self.createRawTable(True))+\
      "\n"+\
      "<p>Legend</p>\n"+\
      "\n"+\
      self.createHtmlTableLegend(True)+\
      "\n"+\
      "<p><b>Trilinos Test/Example Package Dependencies</b></p>\n"+\
      "\n"+\
      self.createHtmlFromTable(self.createRawTable(False))+\
      "\n"+\
      "<p>Legend</p>\n"+\
      "\n"+\
      self.createHtmlTableLegend(False)

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
      "<p><center><b><huge>Trilinos Package Dependencies</huge></b></center></p>\n"\
      "\n"+\
      self.createFullHtmlForTables()+\
      "\n"+\
      "</body>\n"+\
      "\n"+\
      "</html>\n"

    return htmlText


  def writeFullHtmlPage(self, htmlFileName=defaultTrilinosDepsHtmlOutFile):
    htmlFile = open(htmlFileName, 'w')
    htmlFile.write(self.createFullHtmlPage())
    htmlFile.close()


#
# Read in the dependencies from XML
#


def getDependenciesByType(packageEle, typeName):
  packageDepsStr = packageEle.getElementsByTagName(typeName)[0].getAttribute('value');
  if len(packageDepsStr) == 0:
    return []
  return packageDepsStr.split(',')


def getTrilinosDependenciesFromXmlFile(xmlFile=defaultTrilinosDepsXmlInFile):
  #print "xmlFile =", xmlFile
  packageDepXmlDom = xml.dom.minidom.parse(xmlFile)
  trilinosDependencies = TrilinosDependencies()
  for ele in packageDepXmlDom.childNodes[0].childNodes:
    if ele.nodeType == ele.ELEMENT_NODE:
      packageName = ele.getAttribute('name')
      #print "\npackageName =", packageName
      packageDeps = PackageDependencies(packageName,
        getDependenciesByType(ele, "LIB_REQUIRED_DEP_PACKAGES"),
        getDependenciesByType(ele, "LIB_OPTIONAL_DEP_PACKAGES"),
        getDependenciesByType(ele, "TEST_REQUIRED_DEP_PACKAGES"),
        getDependenciesByType(ele, "TEST_OPTIONAL_DEP_PACKAGES")
        )
      #print "\npackageDeps =", str(packageDeps)
      trilinosDependencies.addPackageDependencies(packageDeps)
  return trilinosDependencies
