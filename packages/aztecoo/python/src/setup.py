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

# Certain directory paths are needed by setup.py.  pakDir is the path for the
# epetra package directory, and srcDir is the path for the python source directory
buildDir   = makeInfo.get("top_builddir","")
pakDir     = makeInfo.get("top_srcdir","")
srcDir     = makeInfo.get("srcdir","")

# Define the teuchos include path, library directory and library name
TeuchosSrcDir   = os.path.join(pakDir, "../teuchos", "src")
TeuchosBuildDir = os.path.join(buildDir, "../teuchos", "src")
TeuchosLibDir   = os.path.join("..", "../../teuchos", "src")
TeuchosLib      = "teuchos"

# Define the epetra include path, library directory and library name
EpetraSrcDir   = os.path.join(pakDir, "../epetra", "src")
EpetraBuildDir = os.path.join(buildDir, "../epetra", "src")
EpetraLibDir   = os.path.join("..", "../../epetra", "src")
EpetraLib      = "epetra"
PyEpetraDir   = os.path.join(pakDir, "../epetra/python", "src")


# Define the AztecOO include path, library directory and library name
AztecOOSrcDir   = os.path.join(pakDir, "src")
AztecOOBuildDir = os.path.join(buildDir, "src")
AztecOOLibDir   = os.path.join("..", "..", "src")
AztecOOLib      = "aztecoo"

# Standard libraries.  This is currently a hack.  The library "stdc++" is added
# to the standard library list for a case where we know it needs it.
stdLibs = [ ]
sysName = os.uname()[0]
if sysName == "Linux":
    stdLibs.append("stdc++")

# Create the extra arguments list and complete the standard libraries list.  This
# is accomplished by looping over the arguments in LDFLAGS, FLIBS and LIBS and
# adding them to the appropriate list.
extraArgs = []
libs = makeInfo.get("LDFLAGS"    ,"").split() + \
       makeInfo.get("BLAS_LIBS"  ,"").split() + \
       makeInfo.get("LAPACK_LIBS","").split() + \
       makeInfo.get("FLIBS"      ,"").split() + \
       makeInfo.get("LIBS"       ,"").split()
for lib in libs:
    if lib[:2] == "-l":
        stdLibs.append(lib[2:])
    else:
        extraArgs.append(lib)

# Define the strings that refer to the required source files.
wrapAztecOO          = "AztecOO_wrap.cpp"

# _AztecOO  extension module
_AztecOO = Extension("PyTrilinos._AztecOO",
                    [wrapAztecOO],
                    define_macros=[('HAVE_CONFIG_H', '1')],
                    include_dirs    = [EpetraSrcDir, EpetraBuildDir, 
                                       PyEpetraDir,
                                       TeuchosSrcDir, TeuchosBuildDir, 
                                       AztecOOSrcDir, AztecOOBuildDir, srcDir],
                    library_dirs    = [AztecOOLibDir, TeuchosLibDir, 
                                       EpetraLibDir],
                    libraries       = [AztecOOLib, TeuchosLib, 
                                       EpetraLib] 
                                       + stdLibs,
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
      ext_modules  = [ _AztecOO  ],
      )
