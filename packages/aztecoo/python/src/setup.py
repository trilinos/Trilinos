#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#              PyTrilinos.Amesos : Python Interface to Amesos 
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

# System imports
from   distutils.core import *
from   distutils      import sysconfig
import commands
import os
import string
import sys

# Any information that needs to be transferred from the autotooled Makefile is
# written to file setup.txt using python syntax to define a dictionary.  The
# keys of this 'makeInfo' dictionary are variable names and the corresponding
# values represent the data that will be needed by this setup.py script.
try:
    f = open("setup.txt")
    makeInfo = f.readlines()
    f.close()
    makeInfo = eval(string.join(makeInfo))
except IOError:
    makeInfo = { }

print makeInfo

# Certain directory paths are needed by setup.py.  pakDir is the path for the
# epetra package directory, and srcDir is the path for the python source directory
buildDir   = makeInfo.get("top_builddir","")
pakDir     = makeInfo.get("top_srcdir","")
srcDir     = makeInfo.get("srcdir","")

# Define the epetra include path, library directory and library name
PyEpetraDir   = os.path.join(pakDir, "../epetra/python", "src")

# setup standard information for includes, libraries, and extra agrs
stdIncludes    = [srcDir, PyEpetraDir];
stdLibs        = []
stdLibraryLibs = []
extraArgs      = []

# Create the extra arguments list and complete the standard libraries list.  This
# is accomplished by looping over the arguments in LDFLAGS, FLIBS and LIBS and
# adding them to the appropriate list.
am_libs     = makeInfo.get("LDFLAGS"      ,"").split() + \
              makeInfo.get("AZTECOO_LIBS" ,"").split() + \
              makeInfo.get("BLAS_LIBS"    ,"").split() + \
              makeInfo.get("LAPACK_LIBS"  ,"").split() + \
              makeInfo.get("FLIBS"        ,"").split() + \
              makeInfo.get("LIBS"         ,"").split()
am_includes = makeInfo.get("AZTECOO_INCLUDES" ,"").split()

for lib in am_libs:
    if lib[:2] == "-l":
        stdLibs.append(lib[2:])
    elif lib[:2] == "-L":
        stdLibraryLibs.append(lib[2:])
    else:
        extraArgs.append(lib)

for include in am_includes:
    if include[:2] == "-I":
        stdIncludes.append(include[2:])
    else:
        extraArgs.append(include)

# hack to fix linking under linux
sysName = os.uname()[0]
if sysName == "Linux":
    extraArgs.append("-lstdc++")

# Define the strings that refer to the required source files.
wrapAztecOO          = "AztecOO_wrap.cpp"

# _AztecOO  extension module
_AztecOO = Extension("PyTrilinos._AztecOO",
                    [wrapAztecOO],
                    define_macros=[('HAVE_CONFIG_H', '1')],
                    include_dirs    = stdIncludes,
                    library_dirs    = stdLibraryLibs,
                    libraries       = stdLibs,
                    extra_link_args = extraArgs
                    )

# PyTrilinos.AztecOO setup
setup(name         = "PyTrilinos.AztecOO",
      version      = "1.0",
      description  = "Python Interface to Trilinos Package AztecOO",
      author       = "Marzio Sala",
      author_email = "msala@sandia.gov",
      package_dir  = {"PyTrilinos" : "."},
      packages     = ["PyTrilinos"],
      ext_modules  = [_AztecOO],
      )
