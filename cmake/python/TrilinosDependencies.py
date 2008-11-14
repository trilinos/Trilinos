#!/bin/env python


import xml.dom.minidom


from GeneralScriptSupport import *


class PackageDependencies:
    
  packageName = None
  libRequiredDepPackages = None
  libOptionalDepPackages = None
  testRequiredDepPackages = None
  testOptionalDepPackages = None

  def __init__(self, packageName_in,
    libRequiredDepPackages_in, libOptionalDepPackages_in,
    testRequiredDepPackages_in, testOptionalDepPackages_in
    ):
    self.packageName = packageName_in
    self.libRequiredDepPackages = libRequiredDepPackages_in
    self.libOptionalDepPackages = libOptionalDepPackages_in
    self.testRequiredDepPackages = testRequiredDepPackages_in
    self.testOptionalDepPackages = testOptionalDepPackages_in

  def __str__(self):
    return "{\n"+\
      "  packageName="+self.packageName+",\n"+\
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
    self.__packagesDirToList.update(
      { packageName : len(self.__packagesList)-1 } )

  def __str__(self):
    strRep = ""
    for packageDep in self.__packagesList:
      strRep += str(packageDep)
    return strRep


def getDependenciesByType(packageEle, typeName):
  packageDepsStr = packageEle.getElementsByTagName(typeName)[0].getAttribute('value');
  if len(packageDepsStr) == 0:
    return []
  return packageDepsStr.split(',')


defaultXmlFile = getScriptsDir()+"/data/TrilinosPackageDependencies.xml"


def getTrilinosDependenciesFromXmlFile(xmlFile=defaultXmlFile):
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
