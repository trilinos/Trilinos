#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#              PyTrilinos.Epetra: Python Interface to Epetra
#                   Copyright (2005) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
# ************************************************************************
# @HEADER

# This is not a standard setup.py script that imports distutils, defines
# Extension objects and then calls the setup() function, because the other
# Trilinos packages should have taken care of these steps already.  It does,
# however, emulate setup.py scripts by accepting 'build', 'install', 'clean' and
# 'uninstall' command-line commands.  The first purpose of this script is to
# provide a __init__.py file that specifies all of the PyTrilinos modules.  The
# second purpose is to create shared versions of static Trilinos libraries and
# re-link the Trilinos python extension modules to them.

# Module names: define the names of the modules that could be configured to be a
# part of the PyTrilinos package.  These should be put in dependency order.
modules = [
           "Teuchos",
           #"Thyra",
           "Epetra",
           "TriUtils",
           "EpetraExt",
           "AztecOO",
           "Galeri",
           "Amesos",
           "IFPACK",
           #"Anasazi",
           "ML",
           #"NOX",
           #"LOCA",
           "New_Package",
           ]

# System imports
from   distutils.core import *
from   distutils      import sysconfig
from   distutils.util import get_platform
import os
import sys

# Trilinos import
TRILINOS_HOME_DIR = os.path.normpath(open("TRILINOS_HOME_DIR").read()[:-1])
sys.path.insert(0,os.path.join(TRILINOS_HOME_DIR,"commonTools","buildTools"))
from MakefileVariables import *

# Local import
import SharedUtils

# Build the __init__.py file
def buildInitFile(pyTrilinosModules,filename,depfile):
    build = False
    if os.path.isfile(filename):
        if os.path.getmtime(filename) < os.path.getmtime(depfile):
            build = True
    else:
        build = True
    if build:
        print "creating", filename
        print "\nEnabled modules:"
        for module in pyTrilinosModules:
            print "   ", module
        print
        open(filename,"w").write("__all__ = %s\n" % str(pyTrilinosModules))

# Main script
if __name__ == "__main__":

    # Initialization
    initFileName = "__init__.py"

    # Command-line arguments
    command = sys.argv[1]
    if command not in ("build","install","clean","uninstall"):
        raise RuntimeError, "Command '%s' not supported" % command

    # Determine what packages are enabled
    trilinosExportFile = os.path.join("..","..","..","Makefile.export.trilinos")
    trilinosExport     = processMakefile(trilinosExportFile)
    enabledModules     = [module for module in modules
                          if trilinosExport.has_key("ENABLE_" + module.upper())]

    # Determine the installation prefix
    makeMacros = processMakefile("Makefile")
    prefix     = makeMacros["PYTHON_PREFIX"]
    pyVersion  = "python%d.%d" % sys.version_info[:2]
    installDir = os.path.join(prefix, "lib", pyVersion, "site-packages",
                              "PyTrilinos")

    # Create the shared library builders
    builders = [ ]
    for module in enabledModules:
        builders.append(SharedUtils.SharedTrilinosBuilder(module))

    # Build command
    if command == "build":
        # Convert package libraries to shared
        for builder in builders:
            builder.buildShared()
        # Linke extension modules to dynamic libraries
        for builder in builders:
            builder.reLinkExtension()
        # Build the init file
        buildInitFile(enabledModules,initFileName,trilinosExportFile)

    # Install command
    elif command == "install":
        # Install the shared libraries and extension modules
        for builder in builders:
            builder.install()
        # Build the init file, if needed
        buildInitFile(enabledModules,initFileName,trilinosExportFile)
        # Install
        print "installing", initFileName, "->", installDir
        try:
            open(os.path.join(installDir,initFileName),"w").write(open(initFileName,"r").read())
        except IOError:
            print "\n***Error*** Cannot write to\n%s\n" % installDir
            sys.exit(1)

    # Clean command
    elif command == "clean":
        # Remove any dynamic libraries
        for builder in builders:
            builder.clean()
        # Remove the __init__.py file
        if os.path.isfile(initFileName):
            print "removing", initFileName
            os.remove(initFileName)

    # Uninstall command
    elif command == "uninstall":
        # Uninstall the dynamic libraries
        for builder in builders:
            builder.uninstall()
        # Remove the PyTrilinos package
        print "\nUninstalling PyTrilinos package"
        if os.path.isdir(installDir):
            SharedUtils.runCommand("rm -rf " + installDir)
        else:
            print "nothing needs to be done for", installDir
