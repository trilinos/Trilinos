# @HEADER
# ************************************************************************
#
#             PyTrilinos.AztecOO: Python Interface to AztecOO
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
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
# ************************************************************************
# @HEADER

"""Set the python search path to include the library build directory
created by the python distutils module."""

# System includes
import os
import sys
from   distutils.util import get_platform

# Obtain the current directory name
myDir,myName = os.path.split(__file__)

# Construct the the build library directory name
libDir = "lib.%s-%s" % (get_platform(), sys.version[0:3])

# Get the full path to the build directories
fullPath   = os.path.normpath(os.path.join(myDir, "..", "..", "..", "..", 
  "..", "SERIAL", "packages", "ml", "python", "src", "build", libDir, "PyTrilinos"))
epetraPath = os.path.normpath(os.path.join(myDir, "..", "..", "..", "..",
  "..", "SERIAL", "packages", "epetra", "python", "src", "build", libDir, "PyTrilinos"))
epetraextPath = os.path.normpath(os.path.join(myDir, "..", "..", "..", "..",
  "..", "SERIAL", "packages", "epetraext", "python", "src", "build", libDir, "PyTrilinos"))
aztecooPath = os.path.normpath(os.path.join(myDir, "..", "..", "..", "..",
  "..", "SERIAL", "packages", "aztecoo", "python", "src", "build", libDir, "PyTrilinos"))
galeriPath = os.path.normpath(os.path.join(myDir, "..", "..", "..", "..",
  "..", "SERIAL", "packages", "galeri", "python", "src", "build", libDir, "PyTrilinos"))
ifpackPath = os.path.normpath(os.path.join(myDir, "..", "..", "..", "..",
  "..", "SERIAL", "packages", "ifpack", "python", "src", "build", libDir, "PyTrilinos"))
amesosPath = os.path.normpath(os.path.join(myDir, "..", "..", "..", "..",
  "..", "SERIAL", "packages", "amesos", "python", "src", "build", libDir, "PyTrilinos"))


# Insert the full path to the build library directory
# at the beginning of the python search path
if fullPath:
    sys.path.insert(0,fullPath)
if epetraPath:
    sys.path.insert(1,epetraPath)
if epetraextPath:
    sys.path.insert(2,epetraextPath)
if aztecooPath:
    sys.path.insert(3,aztecooPath)
if galeriPath:
    sys.path.insert(4,galeriPath)
if ifpackPath:
    sys.path.insert(5,ifpackPath)
if amesosPath:
    sys.path.insert(6,amesosPath)
