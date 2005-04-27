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

"""Set the python search path to include the library build directory
created by the python distutils module."""

# System includes
import commands
import os
import string
import sys

# Obtain the dictionary of make information
setup_filename = "../src/setup.txt"
try:
    f = open(setup_filename)
    makeInfo = f.readlines()
    f.close()
    makeInfo = eval(string.join(makeInfo))
except IOError:
    print "WARNING: %s not found" % setup_filename
    makeInfo = { }

# Build the command to get the build library name
cmd = "%s %s/../PyTrilinos/pyLocate --build" % (makeInfo.get("PYTHON"    ,""),
                                                makeInfo.get("top_srcdir",""))
(status,output) = commands.getstatusoutput(cmd)
if status != 0:
    raise RuntimeError, "\n\tUNIX command '%s' gives\n\t%s" % (cmd,output)

# Get the path to the build directory
libDir = os.path.join("..", "src", output, "PyTrilinos")

# Insert the library directory name at the beginning of
# the python search path
if libDir:
    sys.path.insert(0,libDir)
