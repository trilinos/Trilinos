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
pakDir = makeInfo.get("top_srcdir","")
srcDir = makeInfo.get("srcdir"    ,"")
CC     = makeInfo.get("CC" )
CXX    = makeInfo.get("CXX")

# Obtain the version from the package version function, using regular
# expressions.  This assumes that the function returns a string constant of the
# form "LOCA Version xxx - mm/dd/yyyy" and extracts the xxx (which does
# not have to be three characters long).  I know LOCA does not have this yet,
# but maybe someday...
versionRE     = re.compile(r"return.*Version\s+(.*)\s+-\s+\d")
versionHeader = os.path.join(pakDir,"src-loca","src","LOCA_Version.C")
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
teuchosLibDir = os.path.join("..", "..", "..", "..", "teuchos", "src")
teuchosLib    = "teuchos"

# # Define the epetra include path, library directory and library name
# epetraInc    = os.path.join(pakDir, "..", "epetra", "src")
# epetraPyInc  = os.path.join(pakDir, "..", "epetra", "python", "src")
# epetraLibDir = os.path.join("..", "..", "..", "..", "epetra", "src")
# epetraLib    = "epetra"

# Define the nox include path, library directory and library name
noxInc      = os.path.join(pakDir, "src")
noxPyInc    = os.path.join(pakDir, "python", "src")
noxLibDir   = os.path.join("..", "..", "..", "src")
noxLib      = "nox"

# Define the loca include path, library directory and library name
locaInc    = os.path.join(pakDir, "src-loca", "src")
locaLibDir = os.path.join("..", "..", "..", "src-loca", "src")
locaLib    = "loca"

# # Define the nox-epetra include path, library directory and library name
# noxEpetraInc    = os.path.join(pakDir, "src-epetra")
# noxEpetraLibDir = os.path.join("..", "..", "src-epetra")
# noxEpetraLib    = "noxepetra"

# Define the nox-lapack include path, library directory and library name
noxLAPACKInc    = os.path.join(pakDir, "src-lapack")
noxLAPACKLibDir = os.path.join("..", "..", "..", "src-lapack")
noxLAPACKLib    = "noxlapack"

# # Define the loca-epetra include path, library directory and library name
# locaEpetraInc    = os.path.join(pakDir, "src-loca", "src-epetra")
# locaEpetraLibDir = os.path.join("..", "..", "..", "src-loca", "src-epetra")
# locaEpetraLib    = "locaepetra"

# Define the loca-lapack include path, library directory and library name
locaLAPACKInc    = os.path.join(pakDir, "src-loca", "src-lapack")
locaLAPACKLibDir = os.path.join("..", "..", "..", "src-loca", "src-lapack")
locaLAPACKLib    = "localapack"

# # Define the aztecoo include path, library directory and library name
# aztecooInc    = os.path.join(pakDir, "..", "aztecoo", "src")
# aztecooLibDir = os.path.join("..", "..", "..", "aztecoo", "src")
# aztecooLib    = "aztecoo"

# # Define the ifpack include path, library directory and library name
# ifpackInc    = os.path.join(pakDir, "..", "ifpack", "src")
# ifpackLibDir = os.path.join("..", "..", "..", "ifpack", "src")
# ifpackLib    = "ifpack"

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

# Compiler and linker
sysconfig.get_config_vars()
config_vars = sysconfig._config_vars;
config_vars['CC']  = CC
config_vars['CXX'] = CXX

# Define the strings that refer to the required source files.
locaTopLevelWrap          = "LOCA_TopLevel_wrap.cpp"
locaContinuationWrap      = "LOCA_Continuation_wrap.cpp"
locaMultiContinuationWrap = "LOCA_MultiContinuation_wrap.cpp"
locaHomotopyWrap          = "LOCA_Homotopy_wrap.cpp"
locaTimeDependentWrap     = "LOCA_TimeDependent_wrap.cpp"
locaBifurcationWrap       = "LOCA_Bifurcation_wrap.cpp"
locaAbstractWrap          = "LOCA_Abstract_wrap.cpp"
#locaEpetraWrap            = "LOCA_Epetra_wrap.cpp"
locaLAPACKWrap            = "LOCA_LAPACK_wrap.cpp"
locaChanWrap              = "LOCA_Chan_wrap.cpp"
epetraVectorHelper       = os.path.join(srcDir, "Epetra_VectorHelper.cpp")
pyInterface              = os.path.join(srcDir, "PyInterface.cpp"        )

# LOCA_TopLevel extension module
LOCA_TopLevel = Extension("PyTrilinos.LOCA._TopLevel",
                          [locaTopLevelWrap],
                          include_dirs       = [locaInc, teuchosInc, noxInc,
                                                noxPyInc, srcDir],
                          library_dirs       = [locaLibDir, teuchosLibDir,
                                                noxLibDir],
                          libraries          = [locaLib, teuchosLib, noxLib] + \
                                                stdLibs,
                          extra_compile_args = extraCompileArgs,
                          extra_link_args    = extraLinkArgs
                          )

# LOCA_Continuation extension module
LOCA_Continuation = Extension("PyTrilinos.LOCA._Continuation",
                              [locaContinuationWrap],
                              include_dirs       = [locaInc, noxInc, srcDir],
                              library_dirs       = [locaLibDir, noxLibDir],
                              libraries          = [locaLib, noxLib] + stdLibs,
                              extra_compile_args = extraCompileArgs,
                              extra_link_args    = extraLinkArgs
                              )

# LOCA_MultiContinuation extension module
LOCA_MultiContinuation = Extension("PyTrilinos.LOCA._MultiContinuation",
                                   [locaMultiContinuationWrap],
                                   include_dirs       = [locaInc, teuchosInc,
                                                         noxInc, srcDir],
                                   library_dirs       = [locaLibDir],
                                   libraries          = [locaLib] + stdLibs,
                                   extra_compile_args = extraCompileArgs,
                                   extra_link_args    = extraLinkArgs
                                   )

# LOCA_Homotopy extension module
LOCA_Homotopy = Extension("PyTrilinos.LOCA._Homotopy",
                          [locaHomotopyWrap],
                          include_dirs       = [locaInc, noxInc, srcDir],
                          library_dirs       = [locaLibDir],
                          libraries          = [locaLib] + stdLibs,
                          extra_compile_args = extraCompileArgs,
                          extra_link_args    = extraLinkArgs
                          )

# LOCA_TimeDependent extension module
LOCA_TimeDependent = Extension("PyTrilinos.LOCA._TimeDependent",
                               [locaTimeDependentWrap],
                               include_dirs       = [locaInc, noxInc, srcDir],
                               library_dirs       = [locaLibDir],
                               libraries          = [locaLib] + stdLibs,
                               extra_compile_args = extraCompileArgs,
                               extra_link_args    = extraLinkArgs
                               )

# LOCA_Bifurcation extension module
LOCA_Bifurcation = Extension("PyTrilinos.LOCA._Bifurcation",
                             [locaBifurcationWrap],
                             include_dirs       = [locaInc, teuchosInc, noxInc,
                                                   srcDir],
                             library_dirs       = [locaLibDir, noxLibDir],
                             libraries          = [locaLib, noxLib] + stdLibs,
                             extra_compile_args = extraCompileArgs,
                             extra_link_args    = extraLinkArgs
                             )

# LOCA_Abstract extension module
LOCA_Abstract = Extension("PyTrilinos.LOCA._Abstract",
                          [locaAbstractWrap],
                          include_dirs       = [locaInc, teuchosInc, noxInc,
                                                srcDir],
                          library_dirs       = [locaLibDir],
                          libraries          = [locaLib] + stdLibs,
                          extra_compile_args = extraCompileArgs,
                          extra_link_args    = extraLinkArgs
                          )

# # LOCA_Epetra extension module
# LOCA_Epetra = Extension("PyTrilinos.LOCA._Epetra",
#                         #[locaEpetraWrap, callback, epetraVectorHelper,
#                         # numPyArray, numPyWrapper, pyInterface],
#                         [locaEpetraWrap],
#                         include_dirs       = [locaEpetraInc, locaInc, epetraInc,
#                                               epetraPyInc, srcDir],
#                         library_dirs       = [locaEpetraLibDir, locaLibDir,
#                                               aztecooLibDir, ifpackLibDir,
#                                               epetraLibDir],
#                         libraries          = [locaEpetraLib, locaLib, aztecooLib,
#                                               ifpackLib, epetraLib] + stdLibs,
#                         extra_compile_args = extraCompileArgs,
#                         extra_link_args    = extraLinkArgs
#                         )

# LOCA_LAPACK extension module
LOCA_LAPACK = Extension("PyTrilinos.LOCA._LAPACK",
                        [locaLAPACKWrap],
                        include_dirs       = [locaLAPACKInc, noxLAPACKInc,
                                              locaInc, noxInc, teuchosInc,
                                              srcDir],
                        library_dirs       = [locaLAPACKLibDir, noxLAPACKLibDir,
                                              locaLibDir, noxLibDir],
                        libraries          = [locaLAPACKLib, noxLAPACKLib,
                                              locaLib, noxLib] + stdLibs,
                        extra_compile_args = extraCompileArgs,
                        extra_link_args    = extraLinkArgs
                        )

# LOCA_Chan extension module
LOCA_Chan = Extension("PyTrilinos.LOCA._Chan",
                      [locaChanWrap],
                      include_dirs       = [locaInc, srcDir],
                      library_dirs       = [locaLibDir],
                      libraries          = [locaLib] + stdLibs,
                      extra_compile_args = extraCompileArgs,
                      extra_link_args    = extraLinkArgs
                      )

# Build the list of extension modules to wrap
ext_modules = [ LOCA_TopLevel,
                LOCA_Continuation,
                LOCA_MultiContinuation,
                LOCA_Homotopy,
                LOCA_TimeDependent,
                LOCA_Bifurcation,
                LOCA_Abstract           ]
#if makeInfo.get("LOCA_EPETRA",""):
#    ext_modules.append(LOCA_Epetra)
if makeInfo.get("LOCA_LAPACK"):
    ext_modules.append(LOCA_LAPACK)
    #ext_modules.append(LOCA_Chan  )

# PyTrilinos.LOCA setup.  This is what tells the distutils module how to
# create the package.
setup(name         = "PyTrilinos.LOCA",
      version      = version,
      description  = "Python Interface to Trilinos Package LOCA",
      author       = "Bill Spotz",
      author_email = "wfspotz@sandia.gov",
      package_dir  = {"PyTrilinos.LOCA" : "."},
      packages     = ["PyTrilinos", "PyTrilinos.LOCA"],
      ext_modules  = ext_modules
      )
