#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#         PyTrilinos.New_Package: Python Interface to New_Package
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
import re
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
# new_package directory, and srcDir is the path for the python source directory
pakDir = makeInfo.get("top_srcdir","")
srcDir = makeInfo.get("srcdir"    ,"")
CXX    = makeInfo.get("CXX"       ,"")

# Obtain the version from the package version function, using regular
# expressions.  This assumes that the function returns a string constant of the
# form "PackageName Version xxx - mm/dd/yyyy" and extracts the xxx (which does
# not have to be three characters long).
versionRE     = re.compile(r"return.*Version\s+(.*)\s+-\s+\d")
versionHeader = os.path.join(pakDir,"src","New_Package_Version.h")
try:
    header = open(versionHeader)
    lines  = header.readlines()
    header.close()
    for line in lines:
        match = versionRE.search(line)
        if match:
            version = match.group(1)
except IOError:
    version = "??"

# Define the epetra include path, library directory and library name
epetraInc    = os.path.join(pakDir, "..", "epetra", "src")
epetraLibDir = os.path.join("..", "..", "..", "epetra", "src")
epetraLib    = "epetra"

# Define the new_package include path, library directory and library name
newPackageInc    = [os.path.join(pakDir, "src"),
                    os.path.join("..","..","src")]
newPackageLibDir = os.path.join("..", "..", "src")
newPackageLib    = "new_package"

# Standard libraries.  This is currently a bit of a hack.  The library "stdc++"
# is added to the standard library list for a case where we know it sometimes
# needs it.
stdLibs = [ ]
sysName = os.uname()[0]
if sysName == "Linux":
    stdLibs.append("stdc++")

# Create the extra compile and link argument lists and complete the standard
# libraries list.  This is accomplished by looping over the arguments in
# LDFLAGS, BLAS_LIBS, LAPACK_LIBS, FLIBS and LIBS and adding them to the
# appropriate list.
extraCompileArgs = [ ]
if makeInfo.get("NEWP_SWAHILI",""):
    extraCompileArgs.append("-DNEWP_SWAHILI")
extraLinkArgs = [ ]
libs = makeInfo.get("LDFLAGS"    ,"").split() + \
       makeInfo.get("BLAS_LIBS"  ,"").split() + \
       makeInfo.get("LAPACK_LIBS","").split() + \
       makeInfo.get("FLIBS"      ,"").split() + \
       makeInfo.get("LIBS"       ,"").split()
for lib in libs:
    if lib[:2] == "-l":
        stdLibs.append(lib[2:])
    else:
        extraLinkArgs.append(lib)

# Define the strings that refer to the required source files.
newPackageWrap = "New_Package_wrap.cpp"

# Compiler and linker
sysconfig.get_config_vars()
config_vars = sysconfig._config_vars;
config_vars['CC' ] = CXX
config_vars['CXX'] = CXX

# _New_Package extension module
_New_Package = Extension("PyTrilinos._New_Package",
                         [newPackageWrap],
                         include_dirs       = newPackageInc + [epetraInc, srcDir],
                         library_dirs       = [newPackageLibDir, epetraLibDir],
                         libraries          = [newPackageLib, epetraLib] + stdLibs,
                         extra_compile_args = extraCompileArgs,
                         extra_link_args    = extraLinkArgs
                         )

# PyTrilinos.New_Package setup.  This is what tells the distutils module how to
# create the package.
setup(name         = "PyTrilinos.New_Package",
      version      = version,
      description  = "Python Interface to Trilinos Package New_Package",
      author       = "Bill Spotz",
      author_email = "wfspotz@sandia.gov",
      package_dir  = {"PyTrilinos" : "."},
      packages     = ["PyTrilinos"],
      ext_modules  = [ _New_Package ]
      )
