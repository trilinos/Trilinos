#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#                 PyTrilinos.NOX: Python Interface to NOX
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
# nox directory, and srcDir is the path for the python source directory
pakDir   = makeInfo.get("top_srcdir",  "")
srcDir   = makeInfo.get("srcdir",      "")
buildDir = makeInfo.get("top_builddir","")

# Obtain the version from the package version function, using regular
# expressions.  This assumes that the function returns a string constant of the
# form "NOX Version xxx - mm/dd/yyyy" and extracts the xxx (which does
# not have to be three characters long).
versionRE     = re.compile(r"return.*Version\s+(.*)\s+-\s+\d")
versionHeader = os.path.join(pakDir,"src","NOX_Version.C")
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

# Define the teuchos include path, library directory and library name
teuchosInc    = os.path.join(pakDir, "..", "teuchos", "src")
teuchosLibDir = os.path.join("..", "..", "..", "teuchos", "src")
teuchosLib    = "teuchos"

# Define the pytrilinos include path, library directory and library name
pytrilinosInc    = os.path.join(pakDir,   "..", "PyTrilinos", "src")
pytrilinosLibDir = os.path.join(buildDir, "..", "PyTrilinos", "src")
pytrilinosLib    = "pytrilinos"

# Define the epetra include path, library directory and library name
epetraInc    = os.path.join(pakDir,   "..", "epetra", "src")
epetraPyInc  = os.path.join(pakDir,   "..", "epetra", "python", "src")
epetraLibDir = os.path.join(buildDir, "..", "epetra", "src")
epetraLib    = "epetra"

# Define the nox include path, library directory and library name
noxInc    = os.path.join(pakDir,   "src")
noxLibDir = os.path.join(buildDir, "src")
noxLib    = "nox"

# Define the nox-epetra include path, library directory and library name
noxEpetraInc    = os.path.join(pakDir, "src-epetra")
noxEpetraLibDir = os.path.join("..", "..", "src-epetra")
noxEpetraLib    = "noxepetra"

# Define the nox-lapack include path, library directory and library name
noxLAPACKInc    = os.path.join(pakDir, "src-lapack")
noxLAPACKLibDir = os.path.join("..", "..", "src-lapack")
noxLAPACKLib    = "noxlapack"

# Define the aztecoo include path, library directory and library name
aztecooInc    = os.path.join(pakDir, "..", "aztecoo", "src")
aztecooLibDir = os.path.join("..", "..", "..", "aztecoo", "src")
aztecooLib    = "aztecoo"

# Define the ifpack include path, library directory and library name
ifpackInc    = os.path.join(pakDir, "..", "ifpack", "src")
ifpackLibDir = os.path.join("..", "..", "..", "ifpack", "src")
ifpackLib    = "ifpack"

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
extraLinkArgs    = [ ]
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
noxTopLevelWrap    = "NOX_TopLevel_wrap.cpp"
noxAbstractWrap    = "NOX_Abstract_wrap.cpp"
noxParameterWrap   = "NOX_Parameter_wrap.cpp"
noxSolverWrap      = "NOX_Solver_wrap.cpp"
noxStatusTestWrap  = "NOX_StatusTest_wrap.cpp"
noxEpetraWrap      = "NOX_Epetra_wrap.cpp"
noxLAPACKWrap      = "NOX_LAPACK_wrap.cpp"
epetraVectorHelper = os.path.join(srcDir, "Epetra_VectorHelper.cpp")
pyInterface        = os.path.join(srcDir, "PyInterface.cpp"        )

# NOX_TopLevel extension module
NOX_TopLevel = Extension("PyTrilinos.NOX._TopLevel",
                         [noxTopLevelWrap],
                         include_dirs       = [noxInc, srcDir],
                         library_dirs       = [noxLibDir],
                         libraries          = [noxLib] + stdLibs,
                         extra_compile_args = extraCompileArgs,
                         extra_link_args    = extraLinkArgs
                         )

# NOX_Abstract extension module
NOX_Abstract = Extension("PyTrilinos.NOX._Abstract",
                         [noxAbstractWrap],
                         include_dirs       = [noxInc, srcDir],
                         library_dirs       = [noxLibDir],
                         libraries          = [noxLib] + stdLibs,
                         extra_compile_args = extraCompileArgs,
                         extra_link_args    = extraLinkArgs
                         )

# NOX_Parameter extension module
NOX_Parameter = Extension("PyTrilinos.NOX._Parameter",
                          [noxParameterWrap],
                          include_dirs       = [noxInc, teuchosInc, srcDir],
                          library_dirs       = [noxLibDir, teuchosLibDir],
                          libraries          = [noxLib, teuchosLib] + stdLibs,
                          extra_compile_args = extraCompileArgs,
                          extra_link_args    = extraLinkArgs
                          )

# NOX_Solver extension module
NOX_Solver = Extension("PyTrilinos.NOX._Solver",
                       [noxSolverWrap],
                       include_dirs       = [noxInc, srcDir],
                       library_dirs       = [noxLibDir],
                       libraries          = [noxLib] + stdLibs,
                       extra_compile_args = extraCompileArgs,
                       extra_link_args    = extraLinkArgs
                       )

# NOX_StatusTest extension module
NOX_StatusTest = Extension("PyTrilinos.NOX._StatusTest",
                           [noxStatusTestWrap],
                           include_dirs       = [noxInc, srcDir],
                           library_dirs       = [noxLibDir],
                           libraries          = [noxLib] + stdLibs,
                           extra_compile_args = extraCompileArgs,
                           extra_link_args    = extraLinkArgs
                           )

# NOX_Epetra extension module
NOX_Epetra = Extension("PyTrilinos.NOX._Epetra",
                       [noxEpetraWrap, epetraVectorHelper, pyInterface],
                       include_dirs       = [noxEpetraInc, noxInc, epetraInc,
                                             epetraPyInc, pytrilinosInc,
                                             srcDir],
                       library_dirs       = [noxEpetraLibDir, noxLibDir,
                                             aztecooLibDir, ifpackLibDir,
                                             teuchosLibDir, pytrilinosLibDir,
                                             epetraLibDir],
                       libraries          = [noxEpetraLib, noxLib, aztecooLib,
                                             ifpackLib, teuchosLib,
                                             pytrilinosLib, epetraLib] + stdLibs,
                       extra_compile_args = extraCompileArgs,
                       extra_link_args    = extraLinkArgs
                       )

# NOX_LAPACK extension module
NOX_LAPACK = Extension("PyTrilinos.NOX._LAPACK",
                       [noxLAPACKWrap],
                       include_dirs       = [teuchosInc, noxLAPACKInc, noxInc,
                                             srcDir],
                       library_dirs       = [teuchosLibDir, noxLAPACKLibDir,
                                             noxLibDir],
                       libraries          = [teuchosLib, noxLAPACKLib, noxLib] + \
                                            stdLibs,
                       extra_compile_args = extraCompileArgs,
                       extra_link_args    = extraLinkArgs
                       )

# Build the list of extension modules to wrap
ext_modules = [ NOX_TopLevel, NOX_Abstract, NOX_Parameter, NOX_Solver,
                NOX_StatusTest ]
if makeInfo.get("NOX_EPETRA",""):
    ext_modules.append(NOX_Epetra)
if makeInfo.get("NOX_LAPACK",""):
    ext_modules.append(NOX_LAPACK)

# PyTrilinos.NOX setup.  This is what tells the distutils module how to
# create the package.
setup(name         = "PyTrilinos.NOX",
      version      = version,
      description  = "Python Interface to Trilinos Package NOX",
      author       = "Bill Spotz",
      author_email = "wfspotz@sandia.gov",
      package_dir  = {"PyTrilinos.NOX" : "."},
      packages     = ["PyTrilinos", "PyTrilinos.NOX"],
      ext_modules  = ext_modules
      )
