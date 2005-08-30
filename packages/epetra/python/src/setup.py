#! /usr/bin/env python

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

# System imports
from   distutils.core import *
from   distutils      import sysconfig
import commands
import os
import re
import string
import sys

# Trilinos import
TRILINOS_HOME_DIR = os.path.normpath(open("TRILINOS_HOME_DIR","r").read()[:-1])
sys.path.insert(0,os.path.join(TRILINOS_HOME_DIR,"commonTools","buildTools"))

# Build the makeVars dictionary
makeVars = { }
from MakefileVariables import *
makeVars.update(processFile(os.path.join("..","..","Makefile")))
makeVars.update(processFile(os.path.join("..","..","Makefile.export")))
makeVars.update(processFile(os.path.join("..","..","Makefile.export.epetra")))
makeVars.update(processFile(os.path.join("Makefile")))

# Certain directory paths are needed by setup.py.  top_srcdir is the path for the
# epetra package directory, and srcdir is the path for the python source directory
top_srcdir   = os.path.normpath(makeVars["top_srcdir"  ])
srcdir       = os.path.normpath(makeVars["srcdir"      ])
top_builddir = os.path.normpath(makeVars["top_builddir"])

# Obtain the version from the package version function, using regular
# expressions.  This assumes that the function returns a string constant of the
# form "PackageName Version xxx - mm/dd/yyyy" and extracts the xxx (which does
# not have to be three characters long).
version       = "??"
versionRE     = re.compile(r"return.*Version\s+(.*)\s+-\s+\d")
versionHeader = os.path.join(top_srcdir,"src","Epetra_Version.h")
try:
    for line in open(versionHeader).readlines():
        match = versionRE.search(line)
        if match:
            version = match.group(1)
except IOError:
    pass

# Define the PyTrilinos include path, library directory and library name
pytrilinosInc    = os.path.normpath(os.path.join(top_srcdir,   "..", "PyTrilinos", "src"))
pytrilinosLibDir = os.path.normpath(os.path.join(top_builddir, "..", "PyTrilinos", "src"))
pytrilinosLib    = "pytrilinos"

# Start defining arguments that will be needed by the Extension class
include_dirs    = [pytrilinosInc, srcdir ]
library_dirs    = [pytrilinosLibDir      ]
libraries       = [pytrilinosLib         ]
extra_link_args = [                      ]

# Get the relevant Makefile variable values, split them into lists of strings,
# and add them together to to obtain a big list of options
options = makeVars["EPETRA_INCLUDES"].split() + \
          makeVars["EPETRA_LIBS"    ].split()

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

# The library "stdc++" is added to the libraries argument list for a case where
# we know it needs it.  Has this been fixed by forcing the compiler to be CXX?
#sysName = os.uname()[0]
#if sysName == "Linux":
#    libraries.append("stdc++")

# Remove duplicate entries from the argument lists
uniquifyList(include_dirs   )
uniquifyList(library_dirs   )
uniquifyList(libraries      )
uniquifyList(extra_link_args)

# Define the strings that refer to the required source files.
wrapEpetra         = "Epetra_wrap.cpp"
epetraNumPyVector  = os.path.join(srcdir,"Epetra_NumPyVector.cpp" )

# Compiler and linker
CXX = makeVars["CXX"]
sysconfig.get_config_vars()
config_vars = sysconfig._config_vars;
config_vars["CC" ] = CXX
config_vars["CXX"] = CXX

# _Epetra extension module
_Epetra = Extension("PyTrilinos._Epetra",
                    [wrapEpetra, epetraNumPyVector],
                    define_macros   = [("HAVE_CONFIG_H", "1")],
                    include_dirs    = include_dirs,
                    library_dirs    = library_dirs,
                    libraries       = libraries,
                    extra_link_args = extra_link_args
                    )

# PyTrilinos.Epetra setup
setup(name         = "PyTrilinos.Epetra",
      version      = version,
      description  = "Python Interface to Trilinos Package Epetra",
      author       = "Bill Spotz",
      author_email = "wfspotz@sandia.gov",
      package_dir  = {"PyTrilinos" : "."},
      packages     = ["PyTrilinos"],
      ext_modules  = [ _Epetra  ]
      )
