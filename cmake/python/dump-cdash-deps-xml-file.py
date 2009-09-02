#!/usr/bin/env python


#
# Read in the commandline arguments
#

usageHelp = r"""dump-package-deps-table.py [OPTIONS]

Tool that dumps an XML file that can be read by CTest/CDash to specify the
Trilinos package dependenices in a way that is independent of Trilinos.  In
CTest/CDash terminology, a Trilinos package is a "Subproject".

By default, if you just run:

   $ SOME_DIR/dump-cdash-deps-xml-file.py

then the XML file will get written into the main Trilinos source directory
where it can be checked in on the next checkin.

You can also change what XML input file is used and what XML file is written.
This is maily to facilitate unit testing of this script.

Have fun looking through all of the Trilinos dependencies!

"""


from TrilinosDependencies import defaultTrilinosDepsXmlInFile, \
  defaultTrilinosDepsHtmlOutFile, defaultCDashDepsXmlFile

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--input-xml-deps-file", dest="inputXmlDepsFile", type="string",
  default=defaultTrilinosDepsXmlInFile,
  help="Input XML file giving the Trilinos dependencies "+\
    "(default = "+defaultTrilinosDepsXmlInFile+")." )

clp.add_option(
  "--output-cdash-deps-xml-file", dest="outputCDashDepsXmlFile", type="string",
  default=defaultCDashDepsXmlFile,
  help="Output XML file giving the Trilinos dependices CDash language"+\
    "(default = "+defaultCDashDepsXmlFile+")." )

(options, args) = clp.parse_args()


#
# Execute the commands
#


from TrilinosDependencies import getTrilinosDependenciesFromXmlFile

trilinosDependencies = getTrilinosDependenciesFromXmlFile(
  options.inputXmlDepsFile)

trilinosDependencies.writeCDashXmlDepsFile(
  options.outputCDashDepsXmlFile)
