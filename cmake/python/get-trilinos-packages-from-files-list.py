#!/usr/bin/env python

from TrilinosPackageFilePathUtils import *


#
# Read in the commandline arguments
#

usageHelp = r"""get-trilinos-packages-from-files-list.py --file=<FILES_LIST_FILE>
"""

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--files-list-file", dest="filesListFile", type="string", default=None,
  help="File containing the list of modified files." )

clp.add_option(
  "--deps-xml-file", dest="depsXmlFile", type="string",
  default=defaultTrilinosDepsXmlInFile,
  help="File containing the listing of packages, dir names, dependencies, etc."\
    +" (default = '"+defaultTrilinosDepsXmlInFile+"'" )

(options, args) = clp.parse_args()

if not options.filesListFile:
  raise Exception("Error, the option --files-list-file=FILENAME must be set!")

filesList = readStrFromFile(options.filesListFile).split('\n')

trilinosDependencies = getTrilinosDependenciesFromXmlFile(options.depsXmlFile)

packagesList = getPackagesListFromFilePathsList(trilinosDependencies, filesList, True)

print ';'.join(packagesList)
