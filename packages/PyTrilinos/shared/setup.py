#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#                PyTrilinos: Python Interface to Trilinos
#                   Copyright (2007) Sandia Corporation
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

# The purpose of this script is to create shared versions of static Trilinos
# libraries.

# Package names: define the names of the Trilinos pacakages that should be
# converted to shared.  These should be in build order.
packages = [
            "teuchos",
            "rtop",
            "kokkos",
            "epetra",
            "triutils",
            "tpetra",
            "epetraext",
            "thyra",
            "isorropia",
            "aztecoo",
            "galeri",
            "amesos",
            "ifpack",
            "komplex",
            "claps",
            "belos",
            "anasazi",
            "pliris",
            "ml",
            "moertel",
            "stratimikos",
            "meros",
            "sacado",
            "nox",
            "loca",
            "rythmos",
            "moocho",
            ]

# System imports
import os
import sys

# Set the python path to find the PyTrilinos utilities
sys.path.append(os.path.join("..", "util"))

# Local imports
from   MakefileVariables import *
import SharedUtils

# Main logic
def main(command, destdir):

    # Command-line arguments
    if command not in ("build","install","clean","uninstall"):
        raise RuntimeError, "Command '%s' not supported" % command

    # Extract the Makefile variables
    print "\nExtracting Makefile variables ...",
    sys.stdout.flush()
    makeMacros = processMakefile("Makefile")
    print "done"

    # Determine what packages are enabled
    enabledPackages = [package for package in packages
                       if makeMacros.get("ENABLE_" + package.upper(),
                                         "false") == "true"]

    # Determine the installation information
    mkdir      = makeMacros["mkdir_p"]
    libDir     = makeMacros["libdir"]
    if destdir:
        installDir = os.path.join(destdir,*libDir.split(os.sep))
    else:
        installDir = libDir

    ######################################################
    # Build/clean/install/uninstall the shared libraries #
    ######################################################

    if SharedUtils.buildSharedLibraries():

        # Create the shared library builders
        builders = [ ]
        for package in enabledPackages:
            builders.append(SharedUtils.SharedTrilinosBuilder(package))

        # Special cases
        if makeMacros["ENABLE_THYRA" ] == "true" and \
           makeMacros["ENABLE_EPETRA"] == "true":
            builders.append(SharedUtils.SharedTrilinosBuilder("thyraepetra"))
        if makeMacros["ENABLE_THYRA"    ] == "true" and \
           makeMacros["ENABLE_EPETRAEXT"] == "true":
            builders.append(SharedUtils.SharedTrilinosBuilder("thyraepetraext"))
        if makeMacros["ENABLE_STRATIMIKOS"] == "true" and \
           makeMacros["ENABLE_AMESOS"     ] == "true":
            builders.append(SharedUtils.SharedTrilinosBuilder("stratimikosamesos"))
        if makeMacros["ENABLE_STRATIMIKOS"] == "true" and \
           makeMacros["ENABLE_AZTECOO"    ] == "true":
            builders.append(SharedUtils.SharedTrilinosBuilder("stratimikosaztecoo"))
        if makeMacros["ENABLE_STRATIMIKOS"] == "true" and \
           makeMacros["ENABLE_IFPACK"     ] == "true":
            builders.append(SharedUtils.SharedTrilinosBuilder("stratimikosifpack"))
        if makeMacros["ENABLE_STRATIMIKOS"] == "true" and \
           makeMacros["ENABLE_ML"         ] == "true":
            builders.append(SharedUtils.SharedTrilinosBuilder("stratimikosml"))
        if makeMacros["ENABLE_NOX_EPETRA"] == "true":
            builders.append(SharedUtils.SharedTrilinosBuilder("noxepetra"))
        if makeMacros["ENABLE_LOCA"  ] == "true" and \
           makeMacros["ENABLE_EPETRA"] == "true":
            builders.append(SharedUtils.SharedTrilinosBuilder("locaepetra"))
        if makeMacros["ENABLE_MOOCHO"] == "true" and \
           makeMacros["ENABLE_THYRA" ] == "true":
            builders.append(SharedUtils.SharedTrilinosBuilder("moochothyra"))

        # Add the PyTrilinos shared library builder to the end of the list
        if makeMacros["ENABLE_PYTRILINOS"] == "true":
            builders.append(SharedUtils.SharedPyTrilinosBuilder())

        # Build command
        if command == "build":
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
            SharedUtils.runCommand(" ".join([mkdir, installDir]))
            # Install the shared libraries
            for builder in builders:
                builder.install(installDir)

        # Uninstall command
        if command == "uninstall":
            # Uninstall the dynamic libraries
            for builder in builders:
                builder.uninstall()

# Main script
if __name__ == "__main__":
    command = ""
    destdir = ""
    for argument in sys.argv:
        if argument in ("build", "clean", "install", "uninstall"):
            command = argument
        if argument.startswith( "--destdir" ):
            destdir = argument.split( "=" )[1]
    main(command, destdir)
