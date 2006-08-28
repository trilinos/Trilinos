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

# This is not a standard setup.py script, in that it does not depend on the
# distutils module to build python extension modules, because the other Trilinos
# packages should have taken care of that already.  It does, however, emulate
# setup.py scripts by accepting 'build', 'install' and 'clean' command-line
# commands and the '--prefix=...' option for 'install'.  The first purpose of
# this script is to provide a __init__.py file that specifies all of the
# PyTrilinos modules.

# Module names: define the names of the modules that could be configured to be a
# part of the PyTrilinos package.
modules = ["Amesos",
           "Anasazi",
           "AztecOO",
           "Epetra",
           "EpetraExt",
           "Galeri",
           "IFPACK",
           "ML",
           "New_Package",
           "NOX",
           "Teuchos",
           "Thyra",
           "TriUtils",
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

# Build the __init__.py file
def buildInitFile(pyTrilinosModules,filename,depfile):
    build = False
    if os.path.isfile(filename):
        if os.path.getmtime(filename) < os.path.getmtime(depfile):
            build = True
    else:
        build = True
    if build:
        print "\nEnabled modules:"
        for module in pyTrilinosModules:
            print "   ", module
        print
        print "creating", filename
        open(filename,"w").write("__all__ = %s\n" % str(pyTrilinosModules))

# Main script
if __name__ == "__main__":

    # Initialization
    initFileName = "__init__.py"

    # Command-line arguments
    command = sys.argv[1]
    if command not in ("build","install","clean"):
        raise RuntimeError, "Command '%s' not supported" % command

    # Determine what packages are enabled
    trilinosExportFile = os.path.join("..","..","..","Makefile.export.trilinos")
    trilinosExport     = processMakefile(trilinosExportFile)
    enabledModules     = [module for module in modules
                          if trilinosExport.has_key("ENABLE_" + module.upper())]

    # Determine the build directories to copy from
    prefixDir = os.path.abspath(os.path.join("..", "..", "..", "packages"))
    libDir    = "lib.%s-%s" % (get_platform(), sys.version[:3])
    suffixDir = os.path.join("python", "src", "build", libDir, "PyTrilinos")
    copyDirs  = [os.path.join(prefixDir, module.lower(), suffixDir)
                 for module in enabledModules]

    # Build command
    if command == "build":
        buildInitFile(enabledModules,initFileName,trilinosExportFile)

    # Install command
    elif command == "install":
        # Build the install directory name
        if len(sys.argv) > 2 and sys.argv[2].startswith("--prefix="):
            prefix = sys.argv[2].split("=")[1]
        else:
            prefix = sys.prefix
        installDir = os.path.join(prefix,"lib", "python" + sys.version[:3],
                                  "site-packages", "PyTrilinos")
        # Build the init file, if needed
        buildInitFile(enabledModules,initFileName,trilinosExportFile)
        # Install
        print "installing", initFileName, "->", installDir
        try:
            open(os.path.join(installDir,initFileName),"w").write(open(initFileDir,"r").read())
        except:
            print "\n***Error*** Cannot write to\n%s\n" % installDir
            sys.exit(1)

    # Clean command
    elif command == "clean":
        if os.path.isfile(initFileName):
            print "removing", initFileName
            os.remove(initFileName)
