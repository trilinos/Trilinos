#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#                PyTrilinos: Python Interface to Trilinos
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
# Questions? Contact Bill Spotz (wfspotz@sandia.gov)
#
# ************************************************************************
# @HEADER

# The first purpose of this script is to provide a PyTrilinos/__init__.py file
# that specifies all of the PyTrilinos modules.  The second purpose, in the case
# of MPI builds, is to create shared versions of static Trilinos libraries.  The
# third purpose is the standard setup.py purpose: define the distutils Extension
# objects and call the distutils setup function in order to build the PyTrilinos
# package.

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
           "NOX",
           #"LOCA",
           ]

# System imports
from   distutils.core import *
from   distutils      import sysconfig
from   distutils.util import get_platform
import os
import sys

# Local imports
from   MakefileVariables   import *
from   PyTrilinosExtension import *
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
        print "\ncreating", filename
        print "Enabled modules:"
        for module in pyTrilinosModules:
            print "   ", module
        print
        open(filename,"w").write("__all__ = %s\n" % str(pyTrilinosModules))

# Main script
if __name__ == "__main__":

    # Initialization
    initFileName = "PyTrilinos/__init__.py"

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
    print "Extracting Makefile variables ...",
    sys.stdout.flush()
    makeMacros = processMakefile("Makefile")
    print "done"
    prefix     = makeMacros["PYTHON_PREFIX"]
    pyVersion  = "python%d.%d" % sys.version_info[:2]
    installDir = os.path.join(prefix, "lib", pyVersion, "site-packages",
                              "PyTrilinos")

    ######################################################
    # Build/clean/install/uninstall the shared libraries #
    ######################################################

    if SharedUtils.buildSharedLibraries():

        # Create the shared library builders
        builders = [ ]
        for module in enabledModules:
            builders.append(SharedUtils.SharedTrilinosBuilder(module))
        if makeMacros["HAVE_NOX_EPETRA_TRUE"] == "":
            noxEpetraBuilder = SharedUtils.SharedTrilinosBuilder("NoxEpetra")
            builders.append(noxEpetraBuilder)

        # Build command
        if command == "build":
            # Convert package libraries to shared
            for builder in builders:
                builder.buildShared()
            # Link the extension modules to dynamic libraries
            for builder in builders:
                builder.reLinkExtension()

        # Clean command
        elif command == "clean":
            # Remove any dynamic libraries
            for builder in builders:
                builder.clean()

        # Install command
        elif command == "install":
            # Install the shared libraries and extension modules
            for builder in builders:
                builder.install()

        # Uninstall command
        elif command == "uninstall":
            # Uninstall the dynamic libraries
            for builder in builders:
                builder.uninstall()

    #################################################################
    # Build/clean/install/uninstall the PyTrilinos __init__.py file #
    #################################################################

    # Build command
    if command == "build":
        # Build the init file
        buildInitFile(enabledModules,initFileName,trilinosExportFile)

    # Clean command
    elif command == "clean":
        # Remove the __init__.py file
        if os.path.isfile(initFileName):
            print "removing", initFileName
            os.remove(initFileName)

    # Install command
    elif command == "install":
        # Build the init file, if needed
        buildInitFile(enabledModules,initFileName,trilinosExportFile)
        # Install
        print "installing", initFileName, "->", installDir
        try:
            open(os.path.join(installDir,initFileName),"w").write(open(initFileName,"r").read())
        except IOError:
            print "\n***Error*** Cannot write to\n%s\n" % installDir
            sys.exit(1)

    # Uninstall command
    elif command == "uninstall":
        # Remove the PyTrilinos package
        print "\nUninstalling PyTrilinos package"
        if os.path.isdir(installDir):
            SharedUtils.runCommand("rm -rf " + installDir)
        else:
            print "nothing needs to be done for", installDir

    #############################################################
    # Use distutils to build/clean/etc... the extension modules #
    #############################################################

    # Build the list of extension modules
    ext_modules = [ ]
    for module in enabledModules:
        ext_modules.extend(makePyTrilinosExtensions(module))

    # Build the list of package names
    extModNames = [mod.name for mod in ext_modules]
    packages = ["PyTrilinos"]
    for extModName in extModNames:
        if extModName.endswith(".___init__"):
            packages.append(extModName[:-10])

    # Call the distutils setup function.  This defines the PyTrilinos package to
    # distutils and distutils takes over from here.
    setup(name         = "PyTrilinos",
          version      = 4.0,
          description  = "Python interface to Trilinos",
          author       = "Bill Spotz",
          author_email = "wfspotz@sandia.gov",
          packages     = packages,
          ext_modules  = ext_modules
          )
