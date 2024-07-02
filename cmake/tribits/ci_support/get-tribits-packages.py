#!/usr/bin/env python

# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

from FindGeneralScriptSupport import *
from TribitsPackageFilePathUtils import *
from gitdist import addOptionParserChoiceOption


#
# Read in the commandline arguments
#

usageHelp = \
r"""get-tribits-packages.py --deps-xml-file=<DEPS_XML_FILE> \
  --only-top-level-packages=[on|off]

This script returns a comma-separated list of all of the project's top-level
or packages or the full set of packages (i.e. parent and subpackages).
"""

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--deps-xml-file", dest="depsXmlFile", type="string",
  help="File containing TriBITS-generated XML data-structure the listing"+\
    " of packages, dir names, dependencies, etc.")

addOptionParserChoiceOption(
  "--only-top-level-packages", "onlyTopLevelPackagesStr",
  ("on", "off"), 0,
  "If 'on', then only top-level packages will be included.  If 'off', then"+\
  " top-level and subpackages will be included in the list (in order).",
  clp )

(options, args) = clp.parse_args()

if options.onlyTopLevelPackagesStr == "on":
  onlyTopLevelPackages = True
else:
  onlyTopLevelPackages = False

trilinosDependencies = getProjectDependenciesFromXmlFile(options.depsXmlFile)

packagesNamesList = trilinosDependencies.getPackagesNamesList(onlyTopLevelPackages)

print(','.join(packagesNamesList))
