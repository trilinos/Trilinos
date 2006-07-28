#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#                  PyTrilinos.ML: Python Interface to ML
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
import os
import sys

# Trilinos import
TRILINOS_HOME_DIR = os.path.normpath(open("TRILINOS_HOME_DIR").read()[:-1])
sys.path.insert(0,os.path.join(TRILINOS_HOME_DIR,"commonTools","buildTools"))
from MakefileVariables import *

# Build the makeVars dictionary by processing relevant Makefiles
makeVars = { }
makeVars.update(processMakefile("Makefile"))

# Import the variable names and values into the global namespace.  This is
# crucual: every variable name/value pair obtained by processing the specified
# Makefiles above will become actual python variables in the global namespace.
globals().update(makeVars)

# Obtain the package version number string
try:
    version = makeVars["PACKAGE_VERSION"]
except KeyError:
    version = makeVars.get("VERSION","??")

# Initialize arguments that will be needed by the Extension class
include_dirs       = [srcdir]
library_dirs       = [      ]
libraries          = [      ]
extra_link_args    = [      ]
extra_compile_args = CPPFLAGS.split() + CXXFLAGS.split()
uniquifyList(extra_compile_args)

# Get the relevant Makefile export variable values, split them into lists of
# strings, add them together to obtain a big list of option strings, and then
# remove any duplicate entries
options = ML_PYTHON_INCLUDES.split() + \
          ML_PYTHON_LIBS.split()
uniquifyList(options)

# Distribute the individual options to the appropriate Extension class arguments
for option in options:
    if option[:2] == "-I":
        include_dirs.append(option[2:])
    elif option[:2] == "-L":
        library_dirs.append(option[2:])
    elif option[:2] == "-l":
        libraries.append(option[2:])
    else:
        extra_link_args.append(option)

# An additional include directory
include_dirs.append(os.path.join(top_srcdir,"..","teuchos","python","src"))
include_dirs.append(os.path.join(top_srcdir,"..","epetra" ,"python","src"))

# Alter the order of the libraries to ensure proper linking.  This seems to be
# necessary when SuperLU_DIST is enabled for an MPI build and the g2c library is
# required.
g2c = "g2c"
if g2c in libraries:
    libraries.remove(g2c)
    libraries.append(g2c)        

# Define the strings that refer to the required local source files
mlWrap = "ML_wrap.cpp"

# Compiler and linker
sysconfig.get_config_vars()
sysconfig._config_vars["CC" ] = CXX
sysconfig._config_vars["CXX"] = CXX

# _ML extension module
_ML = Extension("PyTrilinos._ML",
                [mlWrap],
                define_macros      = [("HAVE_CONFIG_H", "1"),
                                      ("MLAPI_LC"     , "1")],
                include_dirs       = include_dirs,
                library_dirs       = library_dirs,
                libraries          = libraries,
                extra_compile_args = extra_compile_args,
                extra_link_args    = extra_link_args
                )

# PyTrilinos.ML setup
setup(name         = "PyTrilinos.ML",
      version      = version,
      description  = "Python Interface to Trilinos Package ML",
      author       = "Bill Spotz",
      author_email = "wfspotz@sandia.gov",
      package_dir  = {"PyTrilinos" : "."},
      packages     = ["PyTrilinos"],
      ext_modules  = [ _ML ]
      )
