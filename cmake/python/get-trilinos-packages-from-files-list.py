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

(options, args) = clp.parse_args()

filesList = readStrFromFile(options.filesListFile).split('\n')

trilinosDependencies = getTrilinosDependenciesFromXmlFile(defaultTrilinosDepsXmlInFile)

packagesList = getPackagesListFromFilePathsList(trilinosDependencies, filesList, True)

print ';'.join(packagesList)
