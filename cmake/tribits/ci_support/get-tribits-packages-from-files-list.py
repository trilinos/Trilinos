#!/usr/bin/env python3

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


#
# Read in the commandline arguments
#

usageHelp = \
r"""get-tribits-packages-from-files-list.py --deps-xml-file=<DEPS_XML_FILE> \
    --files-list-file=<FILES_LIST_FILE> [--project-dir=<projectDir>]

This script returns a comma-seprated list of all of the project's TriBITS
packages that must be directly tested for changes in the input list of files.
This may also include the special package name 'ALL_PACKAGES' which means that
at least one changed file (e.g. <projectDir>/CMakeLists.txt) should result in
having to test all of the TriBITS packages in the project.  The logic for
which file changes should trigger testing all packages can be specialized for
the project through the Python module:

  <projectDir>/cmake/ProjectCiFileChangeLogic.py

(if that file exists).

This script is used in continuous integration testing workflows involving
TriBITS projects where only packages impacted by the changes are tested.  For
such a scenario, the list of changed files can come from:

  git diff --name-only <upstream>..<branch-tip>  >  changed-files.txt

where <upstream> (e.g. origin/master) is the commit reference that the local
branch was created from and <branch-tip> is the tip of the topic branch.
"""

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--deps-xml-file", dest="depsXmlFile", type="string",
  help="File containing TriBITS-generated XML data-structure the listing of packages, dir names, dependencies, etc.")

clp.add_option(
  "--files-list-file", dest="filesListFile", type="string", default=None,
  help="File containing the list of modified files relative to project base directory, one file per line." )

clp.add_option(
  "--project-dir", dest="projectDir", type="string", default="",
  help="Base project directory.  Used to access more specialized logic beyond" \
    +" what is known in the <DEPSXMLFILE>.  If empty '', then it will be set" \
    +" automatically if TriBITS is is the standard location w.r.t. the project" \
    +" in relation to this script run from the TriBITS dir.")

(options, args) = clp.parse_args()

if not options.filesListFile:
  raise Exception("Error, the option --files-list-file=FILENAME must be set!")

if not options.projectDir and os.path.isfile(defaultProjectDir+"/ProjectName.cmake"):
  options.projectDir = defaultProjectDir

filesList = readStrFromFile(options.filesListFile).splitlines()

trilinosDependencies = getProjectDependenciesFromXmlFile(options.depsXmlFile)

projectCiFileChangeLogic = getProjectCiFileChangeLogic(options.projectDir)

packagesList = getPackagesListFromFilePathsList(trilinosDependencies, filesList, True,
  projectCiFileChangeLogic)

print(','.join(packagesList))
