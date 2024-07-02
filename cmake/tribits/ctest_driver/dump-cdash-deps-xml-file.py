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

Tool that dumps an XML file that can be read by CTest/CDash to specify the
TriBITS project packages in a way that is independent of TriBITS.  In
CTest/CDash terminology, a TriBITS package is a "Subproject".

By default, if you just run:

   $ <thisDir>/dump-cdash-deps-xml-file.py \
       --input-xml-deps-file=<Project>PackageDependencies.xml \
       --output-cdash-deps-xml-file=CDashSubprojectDependencies.xml

then the XML file CDashSubprojectDependencies.xml will get written.

You can also change what XML input file is used and what XML file is written.
This is mainly to facilitate unit testing of this script.

NOTE: Currently this script is gutted to only list out the parent packages
(subprojects) and their regression email addresses.  It does not currently
list any dependent packages.  Currently CDash does not do anything very useful
with that information anyway and actually makes it really hard to see the
actually subproject-specific test results if there are a lot of upstream
packages (because it lists summaries of all of the dependent packages above it
first).

"""

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--input-xml-deps-file", dest="inputXmlDepsFile", type="string",
  help="Input XML file giving the TriBITS Project Package dependencies.")

clp.add_option(
  "--output-cdash-deps-xml-file", dest="outputCDashDepsXmlFile", type="string",
  help="Output XML file giving the TriBITS project dependencies in XML form" \
  + " that CDash understands.")

(options, args) = clp.parse_args()


#
# Execute the commands
#

import os
import sys

ciSupportDir = os.path.join(
  os.path.dirname(os.path.abspath(__file__)),
  "..", "ci_support" )
sys.path = [ciSupportDir] + sys.path

from TribitsDependencies import getProjectDependenciesFromXmlFile

tribitsDependencies = getProjectDependenciesFromXmlFile(
  options.inputXmlDepsFile)

tribitsDependencies.writeCDashXmlDepsFile(
  options.outputCDashDepsXmlFile)
