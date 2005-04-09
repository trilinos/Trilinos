#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#                  PyTrilinos: Rapid Prototyping Package
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

# Build the python library directory name and library name.  These will be used
# when the linker is invoked to create the python extensions.
# pythonDir = [sysconfig.get_config_var('LIBPL'  )      ]
# pythonLib = [sysconfig.get_config_var('LIBRARY')[3:-2]]

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

# Certain directory paths are needed by setup.py.  pyTDir is the path for the
# PyTrilinos package, pakDir is the path for the Trilinos package directory,
# srcDir is the PyTrilinos source directory path, and noxDir is the PyTrilinos
# NOX source directory path.
pyTDir = makeInfo.get("top_srcdir","")
srcDir = makeInfo.get("srcdir"    ,"")
pakDir = os.path.split(pyTDir      )[0]
noxDir = os.path.join( srcDir,"NOX")

# Define the include paths required by various packages.  Each of these is
# defined as a list of a single or more strings.  Thus they can be added
# together to yield a list of multiple strings.
epetraInc    = [os.path.join(pakDir, "epetra",    "src"             )  ]
epetraExtInc = [os.path.join(pakDir, "epetraext", "src",            ), \
                os.path.join(pakDir, "epetraext", "src", "transform"), \
                os.path.join(pakDir, "epetraext", "src", "coloring" )  ]
noxInc       = [os.path.join(pakDir, "nox",       "src"             )  ]
noxEpetraInc = [os.path.join(pakDir, "nox",       "src-epetra"      )  ]

# Define the library search directories needed to link to various package
# libraries.  Each of these is defined as a list of a single string.  Thus they
# can be added together to yield a list of multiple strings.
epetraLibDir    = [os.path.join("..", "..", "epetra",    "src"       )]
epetraExtLibDir = [os.path.join("..", "..", "epetraext", "src"       )]
noxLibDir       = [os.path.join("..", "..", "nox",       "src"       )]
noxEpetraLibDir = [os.path.join("..", "..", "nox",       "src-epetra")]
aztecLibDir     = [os.path.join("..", "..", "aztecoo",   "src"       )]
ifpackLibDir    = [os.path.join("..", "..", "ifpack",    "src"       )]
teuchosLibDir   = [os.path.join("..", "..", "teuchos",   "src"       )]

# Define the library names for various packages.  Each of these is defined as a
# list of a single string.  Thus they can be added together to yield a list of
# multiple strings.
epetraLib    = ["epetra"   ]
epetraExtLib = ["epetraext"]
noxEpetraLib = ["noxepetra"]
noxLib       = ["nox"      ]
aztecLib     = ["aztecoo"  ]
ifpackLib    = ["ifpack"   ]
teuchosLib   = ["teuchos"   ]

# Standard libraries.  This is currently a hack.  The library "stdc++" is added
# to the standard library list for a case where we know it needs it.
stdLibs = [ ]
sysName = os.uname()[0]
if sysName == "Linux":
    stdLibs.append("stdc++")

# Create the extra arguments list and complete the standard libraries list.  This
# is accomplished by looping over the arguments in FLIB and adding them to the
# appropriate list.
extraArgs = [ ]
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

# Define the strings that refer to the specified source files.  These should be
# combined together as needed in a list as the second argument to the Extension
# constructor.
epetraWrap         = "Epetra_wrap.cxx"
epetraExtWrap      = "EpetraExt_wrap.cxx"
noxAbstractWrap    = os.path.join("NOX", "Abstract_wrap.cxx"      )
noxEpetraWrap      = os.path.join("NOX", "Epetra_wrap.cxx"        )
noxParameterWrap   = os.path.join("NOX", "Parameter_wrap.cxx"     )
noxSolverWrap      = os.path.join("NOX", "Solver_wrap.cxx"        )
noxStatusTestWrap  = os.path.join("NOX", "StatusTest_wrap.cxx"    )
epetraNumPyVector  = os.path.join(srcDir,"Epetra_NumPyVector.cxx" )
epetraVectorHelper = os.path.join(srcDir,"Epetra_VectorHelper.cxx")
numPyArray         = os.path.join(srcDir,"NumPyArray.cxx"         )
numPyWrapper       = os.path.join(srcDir,"NumPyWrapper.cxx"       )
noxCallback        = os.path.join(noxDir,"Callback.cxx"           )
noxPyInterface     = os.path.join(noxDir,"PyInterface.cxx"        )

# Epetra extension module
Epetra = Extension("PyTrilinos._Epetra",
                   [epetraWrap,
                    epetraNumPyVector,
                    epetraVectorHelper,
                    numPyArray,
                    numPyWrapper],
                   include_dirs    = epetraInc + [srcDir],
                   library_dirs    = epetraLibDir + teuchosLibDir,
                   libraries       = epetraLib + teuchosLib + stdLibs,
                   extra_link_args = extraArgs
                   )

# EpetraExt extension module
EpetraExt = Extension("PyTrilinos._EpetraExt",
                      [epetraExtWrap],
                      include_dirs    = epetraInc + epetraExtInc,
                      library_dirs    = epetraLibDir + epetraExtLibDir + teuchosLibDir,
                      libraries       = epetraExtLib + epetraLib + teuchosLib  + stdLibs,
                      extra_link_args = extraArgs
                      )

# NOX_Epetra extension module
NOX_Epetra = Extension("PyTrilinos.NOX._Epetra",
                       [noxEpetraWrap,
                        epetraVectorHelper,
                        numPyArray,
                        numPyWrapper,
                        noxCallback,
                        noxPyInterface],
                       include_dirs    = epetraInc + noxInc + noxEpetraInc + \
                                         [srcDir, noxDir],
                       library_dirs    = epetraLibDir + noxLibDir + noxEpetraLibDir + \
                                         aztecLibDir + ifpackLibDir + teuchosLibDir,
                       libraries       = noxEpetraLib + noxLib  + \
                                         aztecLib + ifpackLib + epetraLib + \
                                         teuchosLib + stdLibs,
                       extra_link_args = extraArgs
                       )

# NOX.Abstract extension module
NOX_Abstract = Extension("PyTrilinos.NOX._Abstract",
                         [noxAbstractWrap],
                         include_dirs    = noxInc,
                         library_dirs    = noxLibDir + teuchosLibDir,
                         libraries       = noxLib + teuchosLib + stdLibs,
                         extra_link_args = extraArgs
                         )

# NOX.Parameter extension module
NOX_Parameter = Extension("PyTrilinos.NOX._Parameter",
                          [noxParameterWrap],
                          include_dirs    = noxInc,
                          library_dirs    = noxLibDir + teuchosLibDir,
                          libraries       = noxLib + teuchosLib + stdLibs,
                          extra_link_args = extraArgs
                          )

# NOX.Solver extension module
NOX_Solver = Extension("PyTrilinos.NOX._Solver",
                       [noxSolverWrap],
                       include_dirs    = noxInc,
                       library_dirs    = noxLibDir + teuchosLibDir,
                       libraries       = noxLib + teuchosLib + stdLibs,
                       extra_link_args = extraArgs
                       )

# NOX.StatusTest extension module
NOX_StatusTest = Extension("PyTrilinos.NOX._StatusTest",
                           [noxStatusTestWrap],
                           include_dirs    = noxInc,
                           library_dirs    = noxLibDir + teuchosLibDir,
                           libraries       = noxLib + teuchosLib + stdLibs,
                           extra_link_args = extraArgs
                           )

# PyTrilinos setup
setup(name         = "PyTrilinos",
      version      = "1.0",
      description  = "Python Trilinos Interface",
      author       = "Bill Spotz",
      author_email = "wfspotz@sandia.gov",
      package_dir  = {"PyTrilinos" : "."},
      packages     = ["PyTrilinos", "PyTrilinos.NOX"],
      ext_modules  = [ Epetra,       EpetraExt,     NOX_Epetra,
                       NOX_Abstract, NOX_Parameter, NOX_Solver,
                       NOX_StatusTest                          ]
      )
