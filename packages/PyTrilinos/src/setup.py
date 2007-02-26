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
           "Anasazi",
           "ML",
           "NOX",
           "LOCA",
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
def buildInitFile(filename,depfile,pyTrilinosModules,
                  pyTrilinosVersion,trilinosVersion):
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
        content = """
__all__ = %s
__version__ = '%s'
def version(): return 'Trilinos version: %s\\nPyTrilinos version: ' + __version__
        """ % (str(pyTrilinosModules), pyTrilinosVersion, trilinosVersion)
        open(filename,"w").write(content)

# Main script
if __name__ == "__main__":

    # Initialization
    initFileName = os.path.join("PyTrilinos", "__init__.py")

    # Command-line arguments
    command = sys.argv[1]
    if command not in ("build","install","clean","uninstall"):
        raise RuntimeError, "Command '%s' not supported" % command

    # Extract the Makefile variables
    print "\nExtracting Makefile variables ...",
    sys.stdout.flush()
    makeMacros = processMakefile("Makefile")
    print "done"

    # Determine what packages are enabled
    enabledModules = [module for module in modules
                      if makeMacros["ENABLE_" + module.upper()] == "true"]

    # Determine the Trilinos version
    trilinosVersion = processMakefile(os.path.join("..","..","..","Makefile"))["PACKAGE_VERSION"]

    # Determine the installation information
    srcdir            = makeMacros["srcdir"]
    prefix            = makeMacros["prefix"]
    pythonPrefix      = makeMacros["PYTHON_PREFIX"]
    pyTrilinosVersion = makeMacros["PACKAGE_VERSION"]
    pyVersion         = "python%d.%d" % sys.version_info[:2]
    install           = makeMacros["INSTALL"]
    mkdir             = makeMacros["mkdir_p"]
    libDir            = os.path.join(prefix, "lib")
    pyTrilinosDir     = os.path.join(pythonPrefix, "lib", pyVersion, "site-packages",
                                     "PyTrilinos")

    ######################################################
    # Build/clean/install/uninstall the shared libraries #
    ######################################################

    if SharedUtils.buildSharedLibraries():

        # Create the shared library builders
        builders = [ ]
        for module in enabledModules:
            builders.append(SharedUtils.SharedTrilinosBuilder(module))
        if makeMacros["ENABLE_NOX_EPETRA"] == "true":
            noxEpetraBuilder = SharedUtils.SharedTrilinosBuilder("NoxEpetra")
            builders.append(noxEpetraBuilder)

        # Build command
        if command in ("build", "install"):
            # Convert package libraries to shared
            for builder in builders:
                builder.buildShared()

        # Clean command
        if command == "clean":
            # Remove any dynamic libraries
            for builder in builders:
                builder.clean()

        # Install command
        if command == "install":
            # Make sure the lib directory exists
            SharedUtils.runCommand(" ".join([mkdir, libDir]))
            # Install the shared libraries and extension modules
            for builder in builders:
                builder.install()

        # Uninstall command
        if command == "uninstall":
            # Uninstall the dynamic libraries
            for builder in builders:
                builder.uninstall()

    #########################################################
    # Build/clean/uninstall the PyTrilinos __init__.py file #
    #########################################################

    # Build command
    if command in ("build", "install"):
        # Build the init file
        buildInitFile(initFileName, "Makefile", enabledModules,
                      pyTrilinosVersion, trilinosVersion)

    # Clean command
    if command == "clean":
        # Remove the __init__.py file
        if os.path.isfile(initFileName):
            print "removing", initFileName
            os.remove(initFileName)

    # Uninstall command
    if command == "uninstall":
        # Remove the PyTrilinos package
        print "\nUninstalling PyTrilinos package"
        if os.path.isdir(pyTrilinosDir):
            SharedUtils.runCommand("rm -rf " + pyTrilinosDir)
        else:
            print "nothing needs to be done for", pyTrilinosDir
        # "uninstall" is not a distutils command, so end here
        sys.exit()

    ###################
    # UserArray patch #
    ###################

    # NumPy version 0.9.8 has a bug in UserArray.  If the user is using this
    # version of NumPy, we need to include our patched version of UserArray.py
    # in the distribution
    if command in ("build", "install"):
        from numpy import __version__ as numpy_version
        if numpy_version == "0.9.8":
            userArraySrc = os.path.join(srcdir,"UserArray.patch")
            userArrayTrg = "UserArrayFix.py"
            if SharedUtils.needsToBeBuilt(userArrayTrg, [userArraySrc]):
                print "copying %s -> %s" % (userArraySrc, userArrayTrg)
                open(userArrayTrg,"w").write(open(userArraySrc,"r").read())

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
          version      = pyTrilinosVersion,
          description  = "Python interface to Trilinos",
          author       = "Bill Spotz",
          author_email = "wfspotz@sandia.gov",
          url          = "http://software.sandia.gov/trilinos/packages/pytrilinos",
          download_url = "http://software.sandia.gov/trilinos/downloads/trilinos-7.0.html",
          packages     = packages,
          ext_modules  = ext_modules
          )
