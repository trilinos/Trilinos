#!/usr/bin/env python

# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

from TribitsPackageFilePathUtils import *


#
# Read in the command-line arguments
#

usageHelp = \
r"""filter-packages-list.py --deps-xml-file=<PROJECT_DEPS_FILE> \
  --input-packages-list=<P1>,<P2>,... --keep-test-test-categories=<T1>,<T2>,...

This script takes in a comma-separated list of TriBITS package names in
--input-packages-list='<P1>,<P2>,...' and then filters out the names for
packages that don't match the set categories listed in
--keep-test-test-categories='<T1>,<T2>,...' (where each package's test test
category is given the input TriBITS-generated project dependencies file
--deps-xml-file=<PROJECT_DEPS_FILE>).  The filtered list of packages is
printed to STDOUT as a comma-separated list.

For example, to keep only the Primary Tested (PT) packages, use:

  filter-packages-list.py --keep-test-test-categories=PT [other args]

To keep both Primary Tested and Secondary Tested packages, use:

  filter-packages-list.py --keep-test-test-categories=PT,ST [other args]

To keep all packages, use:

  filter-packages-list.py --keep-test-test-categories=PT,ST,EX [other args]

(or don't bother running the script).
"""

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--deps-xml-file", dest="depsXmlFile", type="string",
  help="TriBITS generated XML file containing the listing of packages, dir names, dependencies, etc.")

clp.add_option(
  "--input-packages-list", dest="inputPackagesList", type="string", default="",
  help="Comma-seprated List of packages that needs to be filtered (i.e. \"P1,P2,...\")." )

clp.add_option(
  "--keep-test-test-categories", dest="keepTestTestCategories", type="string", default="",
  help="List of package types to keep (i.e. \"PT,ST,EX\"." )

(options, args) = clp.parse_args()

trilinosDependencies = getProjectDependenciesFromXmlFile(options.depsXmlFile)

inputPackagesList = options.inputPackagesList.split(",")
keepTestTestCategoriesList = options.keepTestTestCategories.split(",")
outputPackagesList = \
  trilinosDependencies.filterPackageNameList(inputPackagesList, keepTestTestCategoriesList)
print(','.join(outputPackagesList))
