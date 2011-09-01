#!/usr/bin/env python

from TrilinosPackageFilePathUtils import *


#
# Read in the commandline arguments
#

usageHelp = r"""filter-packages-list --input-packages-list=P1,P2,... --keep-types=T1,T2,...
"""

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--input-packages-list", dest="inputPackagesList", type="string", default="",
  help="List of packages that needs to be filtered (i.e. \"P1,P2,...\")." )

clp.add_option(
  "--keep-types", dest="keepTypes", type="string", default="",
  help="List of types to keep (i.e. \"PS,SS,EX\"." )

clp.add_option(
  "--deps-xml-file", dest="depsXmlFile", type="string", default=defaultTrilinosDepsXmlInFile,
  help="File containing the listing of packages, dir names, dependencies, etc."\
    +" (default = '"+defaultTrilinosDepsXmlInFile+"'" )

(options, args) = clp.parse_args()

trilinosDependencies = getTrilinosDependenciesFromXmlFile(options.depsXmlFile)

inputPackagesList = options.inputPackagesList.split(",")
keepTypesList = options.keepTypes.split(",")
outputPackagesList = trilinosDependencies.filterPackageNameList(inputPackagesList, keepTypesList)
print ','.join(outputPackagesList)
