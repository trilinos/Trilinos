#!/usr/bin/env python

# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


#
# Read in the commandline arguments
#

usageHelp = r"""dump-package-deps-table.py [OPTIONS]

Tool that creates an HTML page showing the TriBITS Project dependencies which
are specified in an input XML file.

By default, if you just run:

   $ SOME_DIR/dump-package-deps-table.py

then the table will get written into the main TriBITS Project source directory
where it can be checked in on the next checkin.

You can also change what XML input file is used and what HTML file is written.
This is maining to facilitate unit testing of this code.

Have fun looking through all of the TriBITS Project dependencies!

"""


from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--input-xml-deps-file", dest="inputXmlDepsFile", type="string",
  help="Input XML file giving the TriBITS Project dependencies.")

clp.add_option(
  "--output-html-deps-file", dest="outputHtmlDepsFile", type="string",
  help="Output HTML file showing the TriBITS Project dependencies")

(options, args) = clp.parse_args()


#
# Execute the commands
#


from TribitsDependencies import getProjectDependenciesFromXmlFile

trilinosDependencies = getProjectDependenciesFromXmlFile(
  options.inputXmlDepsFile)

trilinosDependencies.writeFullHtmlPage(
  options.outputHtmlDepsFile)
