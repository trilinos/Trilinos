#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#                PyTrilinos.LOCA: Python Interface to LOCA
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
makeVars.update(processMakefile(os.path.join("Makefile")))

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
include_dirs    = [srcdir]
library_dirs    = [      ]
libraries       = [      ]
extra_link_args = [      ]

# Get the relevant Makefile export variable values, split them into lists of
# strings, add them together to obtain a big list of option strings, and then
# remove any duplicate entries
options = LOCA_PYTHON_INCLUDES.split() + \
          LOCA_PYTHON_LIBS.split()
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
include_dirs.append(os.path.join(top_srcdir,"..","epetra","python","src"))

# Define the strings that refer to the required local source files
locaTopLevelWrap          = "LOCA_TopLevel_wrap.cpp"
locaAbstractWrap          = "LOCA_Abstract_wrap.cpp"
locaBifurcationWrap       = "LOCA_Bifurcation_wrap.cpp"
locaContinuationWrap      = "LOCA_Continuation_wrap.cpp"
locaHomotopyWrap          = "LOCA_Homotopy_wrap.cpp"
locaMultiContinuationWrap = "LOCA_MultiContinuation_wrap.cpp"
locaTimeDependentWrap     = "LOCA_TimeDependent_wrap.cpp"
locaEpetraWrap            = "LOCA_Epetra_wrap.cpp"
locaLAPACKWrap            = "LOCA_LAPACK_wrap.cpp"

# Compiler and linker
sysconfig.get_config_vars()
sysconfig._config_vars["CC" ] = CXX
sysconfig._config_vars["CXX"] = CXX

# LOCA_TopLevel extension module
LOCA_TopLevel = Extension("PyTrilinos.LOCA_TopLevel",
                          [locaTopLevelWrap],
                          define_macros   = [("HAVE_CONFIG_H", "1")],
                          include_dirs    = include_dirs,
                          library_dirs    = library_dirs,
                          libraries       = libraries,
                          extra_link_args = extra_link_args
                          )

# LOCA_Abstract extension module
LOCA_Abstract = Extension("PyTrilinos.LOCA_Abstract",
                          [locaAbstractWrap],
                          define_macros   = [("HAVE_CONFIG_H", "1")],
                          include_dirs    = include_dirs,
                          library_dirs    = library_dirs,
                          libraries       = libraries,
                          extra_link_args = extra_link_args
                          )

# LOCA_Bifurcation extension module
LOCA_Bifurcation = Extension("PyTrilinos.LOCA_Bifurcation",
                             [locaBifurcationWrap],
                             define_macros   = [("HAVE_CONFIG_H", "1")],
                             include_dirs    = include_dirs,
                             library_dirs    = library_dirs,
                             libraries       = libraries,
                             extra_link_args = extra_link_args
                             )

# LOCA_Continuation extension module
LOCA_Continuation = Extension("PyTrilinos.LOCA_Continuation",
                              [locaContinuationWrap],
                              define_macros   = [("HAVE_CONFIG_H", "1")],
                              include_dirs    = include_dirs,
                              library_dirs    = library_dirs,
                              libraries       = libraries,
                              extra_link_args = extra_link_args
                              )

# LOCA_Homotopy extension module
LOCA_Homotopy = Extension("PyTrilinos.LOCA_Homotopy",
                          [locaHomotopyWrap],
                          define_macros   = [("HAVE_CONFIG_H", "1")],
                          include_dirs    = include_dirs,
                          library_dirs    = library_dirs,
                          libraries       = libraries,
                          extra_link_args = extra_link_args
                          )

# LOCA_MultiContinuation extension module
LOCA_MultiContinuation = Extension("PyTrilinos.LOCA_MultiContinuation",
                                   [locaMultiContinuationWrap],
                                   define_macros   = [("HAVE_CONFIG_H", "1")],
                                   include_dirs    = include_dirs,
                                   library_dirs    = library_dirs,
                                   libraries       = libraries,
                                   extra_link_args = extra_link_args
                                   )

# LOCA_TimeDependent extension module
LOCA_TimeDependent = Extension("PyTrilinos.LOCA_TimeDependent",
                               [locaTimeDependentWrap],
                               define_macros   = [("HAVE_CONFIG_H", "1")],
                               include_dirs    = include_dirs,
                               library_dirs    = library_dirs,
                               libraries       = libraries,
                               extra_link_args = extra_link_args
                               )

# LOCA_Epetra extension module
LOCA_Epetra = Extension("PyTrilinos.LOCA_Epetra",
                        [locaEpetraWrap],
                        define_macros   = [("HAVE_CONFIG_H", "1")],
                        include_dirs    = include_dirs,
                        library_dirs    = library_dirs,
                        libraries       = libraries,
                        extra_link_args = extra_link_args
                        )

# LOCA_LAPACK extension module
LOCA_LAPACK = Extension("PyTrilinos.LOCA_LAPACK",
                        [locaLAPACKWrap],
                        define_macros   = [("HAVE_CONFIG_H", "1")],
                        include_dirs    = include_dirs,
                        library_dirs    = library_dirs,
                        libraries       = libraries,
                        extra_link_args = extra_link_args
                        )

# Build the list of extension modules to wrap
ext_modules = [LOCA_TopLevel, LOCA_Abstract, LOCA_Bifurcation,
               LOCA_Continuation, LOCA_Homotopy, LOCA_MultiContinuation,
               LOCA_TimeDependent]
if LOCA_EPETRA: ext_modules.append(LOCA_Epetra)
if LOCA_LAPACK: ext_modules.append(LOCA_LAPACK)

# PyTrilinos.LOCA setup
setup(name         = "PyTrilinos.LOCA",
      version      = version,
      description  = "Python Interface to Trilinos Package LOCA",
      author       = "Bill Spotz",
      author_email = "wfspotz@sandia.gov",
      package_dir  = {"PyTrilinos.LOCA" : "."},
      packages     = ["PyTrilinos", "PyTrilinos.LOCA"],
      ext_modules  = ext_modules
      )
