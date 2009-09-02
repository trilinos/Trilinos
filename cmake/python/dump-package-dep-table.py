#!/usr/bin/env python


#
# Read in the commandline arguments
#

usageHelp = r"""dump-package-deps-table.py [OPTIONS]

Tool that creates an HTML page showing the Trilinos dependencies which are
specified in an input XML file.

By default, if you just run:

   $ SOME_DIR/dump-package-deps-table.py

then the table will get written into the main Trilinos source directory where
it can be checked in on the next checkin.

You can also change what XML input file is used and what HTML file is written.
This is maining to facilitate unit testing of this code.

Have fun looking through all of the Trilinos dependencies!

"""


from TrilinosDependencies import defaultTrilinosDepsXmlInFile, \
  defaultTrilinosDepsHtmlOutFile

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--input-xml-deps-file", dest="inputXmlDepsFile", type="string",
  default=defaultTrilinosDepsXmlInFile,
  help="Input XML file giving the Trilinos dependencies "+\
    "(default = "+defaultTrilinosDepsXmlInFile+")." )

clp.add_option(
  "--output-html-deps-file", dest="outputHtmlDepsFile", type="string",
  default="TrilinosPackageDependenciesTable.html",
  help="Output HTML file showing the Trilinos dependencies"+\
    "(default = TrilinosPackageDependenciesTable.html)." )

(options, args) = clp.parse_args()


#
# Execute the commands
#


from TrilinosDependencies import getTrilinosDependenciesFromXmlFile

trilinosDependencies = getTrilinosDependenciesFromXmlFile(
  options.inputXmlDepsFile)

trilinosDependencies.writeFullHtmlPage(
  options.outputHtmlDepsFile)
