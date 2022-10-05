#!/usr/bin/env python

# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
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
