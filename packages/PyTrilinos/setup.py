#! /usr/bin/env python

# System imports
from   distutils.core import *
from   distutils      import sysconfig
import commands
import os
import string
import sys

# Build the python library directory name and library name
pythonDir = [sysconfig.get_config_var('LIBPL'  )      ]
pythonLib = [sysconfig.get_config_var('LIBRARY')[3:-2]]

# Read the dictionary from file setup.txt
f = open("setup.txt")
makeInfo = f.readlines()
f.close()
makeInfo = eval(string.join(makeInfo))

# Package, PyTrilinos, source and NOX directories
pyTDir = makeInfo['srcdir']
pakDir = os.path.split(pyTDir      )[0]
srcDir = os.path.join( pyTDir,"src")
noxDir = os.path.join( srcDir,"NOX")

# Include directories
epetraInc    = [os.path.join(pakDir, "epetra",    "src"             )  ]
epetraExtInc = [os.path.join(pakDir, "epetraext", "src", "transform"), \
                os.path.join(pakDir, "epetraext", "src", "coloring" )  ]
noxInc       = [os.path.join(pakDir, "nox",       "src"             )  ]
noxEpetraInc = [os.path.join(pakDir, "nox",       "src-epetra"      )  ]

# Library directories
epetraLibDir    = [os.path.join("..", "epetra",    "src"       )]
epetraExtLibDir = [os.path.join("..", "epetraext", "src"       )]
noxLibDir       = [os.path.join("..", "nox",       "src"       )]
noxEpetraLibDir = [os.path.join("..", "nox",       "src-epetra")]
aztecLibDir     = [os.path.join("..", "aztecoo",   "src"       )]
ifpackLibDir    = [os.path.join("..", "ifpack",    "src"       )]

# Library names
epetraLib    = ["epetra"   ]
epetraExtLib = ["epetraext"]
noxEpetraLib = ["noxepetra"]
noxLib       = ["nox"      ]
aztecLib     = ["aztecoo"  ]
ifpackLib    = ["ifpack"   ]

# LAPACK hack
sysName = os.uname()[0]
if sysName == "Darwin":
    lapackFwk = ["-framework", "veclib"]
    lapackLib = [ ]
else:
    lapackFwk = [ ]
    lapackLib = ["blas", "lapack", "g2c"]

# SWIG Python library
swigStrings = commands.getoutput("swig -python -ldflags").split()
swigDir     = [swigStrings[0][2:]]
swigLib     = [swigStrings[1][2:]]
if sysName == "Linux":
    swigPath = ["-Wl,-rpath"] + swigDir
else:
    swigPath = [ ]

# Standard libraries
if sysName == "Linux":
    stdCXX = ["stdc++"]
else:
    stdCXX = [ ]
stdLibs    = stdCXX + swigLib

# Source files
rawEpetraWrap      = os.path.join(srcDir,"RawEpetra_wrap.cxx"     )
epetraExtWrap      = os.path.join(srcDir,"EpetraExt_wrap.cxx"     )
noxAbstractWrap    = os.path.join(noxDir,"Abstract_wrap.cxx"      )
noxEpetraWrap      = os.path.join(noxDir,"Epetra_wrap.cxx"        )
noxParameterWrap   = os.path.join(noxDir,"Parameter_wrap.cxx"     )
noxSolverWrap      = os.path.join(noxDir,"Solver_wrap.cxx"        )
noxStatusTestWrap  = os.path.join(noxDir,"StatusTest_wrap.cxx"    )
epetraNumPyVector  = os.path.join(srcDir,"Epetra_NumPyVector.cxx" )
epetraVectorHelper = os.path.join(srcDir,"Epetra_VectorHelper.cxx")
numPyArray         = os.path.join(srcDir,"NumPyArray.cxx"         )
numPyWrapper       = os.path.join(srcDir,"NumPyWrapper.cxx"       )
noxCallback        = os.path.join(noxDir,"Callback.cxx"           )
noxPyInterface     = os.path.join(noxDir,"PyInterface.cxx"        )

# Epetra extension module
RawEpetra = Extension("PyTrilinos._RawEpetra",
                      [rawEpetraWrap,
                       epetraNumPyVector,
                       epetraVectorHelper,
                       numPyArray,
                       numPyWrapper],
                      include_dirs    = epetraInc,
                      library_dirs    = epetraLibDir + swigDir + pythonDir,
                      libraries       = stdLibs + epetraLib + lapackLib + pythonLib,
                      extra_link_args = swigPath + lapackFwk
                      )

# EpetraExt extension module
EpetraExt = Extension("PyTrilinos._EpetraExt",
                      [epetraExtWrap],
                      include_dirs    = epetraInc + epetraExtInc,
                      library_dirs    = epetraLibDir + epetraExtLibDir + swigDir + \
                                        pythonDir,
                      libraries       = stdLibs + epetraLib + epetraExtLib + \
                                        lapackLib + pythonLib,
                      extra_link_args = swigPath + lapackFwk
                      )

# NOX_Epetra extension module
NOX_Epetra = Extension("PyTrilinos.NOX._Epetra",
                       [noxEpetraWrap,
                        epetraVectorHelper,
                        numPyArray,
                        numPyWrapper,
                        noxCallback,
                        noxPyInterface],
                       include_dirs    = epetraInc + noxInc + noxEpetraInc,
                       library_dirs    = epetraLibDir + noxLibDir + noxEpetraLibDir + \
                                         swigDir,
                       libraries       = stdLibs + epetraLib + noxLib + noxEpetraLib + \
                                         aztecLib + ifpackLib + lapackLib,
                       extra_link_args = swigPath + lapackFwk
                       )

# NOX.Abstract extension module
NOX_Abstract = Extension("PyTrilinos.NOX._Abstract",
                         [noxAbstractWrap],
                         include_dirs    = noxInc,
                         library_dirs    = noxLibDir + swigDir,
                         libraries       = stdLibs + noxLib,
                         extra_link_args = swigPath
                         )

# NOX.Parameter extension module
NOX_Parameter = Extension("PyTrilinos.NOX._Parameter",
                          [noxParameterWrap],
                          include_dirs    = noxInc,
                          library_dirs    = noxLibDir + swigDir,
                          libraries       = stdLibs + noxLib,
                          extra_link_args = swigPath
                          )

# NOX.Solver extension module
NOX_Solver = Extension("PyTrilinos.NOX._Solver",
                       [noxSolverWrap],
                       include_dirs    = noxInc,
                       library_dirs    = noxLibDir + swigDir,
                       libraries       = stdLibs + noxLib,
                       extra_link_args = swigPath
                       )

# NOX.StatusTest extension module
NOX_StatusTest = Extension("PyTrilinos.NOX._StatusTest",
                           [noxStatusTestWrap],
                           include_dirs    = noxInc,
                           library_dirs    = noxLibDir + swigDir,
                           libraries       = stdLibs + noxLib,
                           extra_link_args = swigPath
                           )

# PyTrilinos setup
setup(name         = "PyTrilinos",
      version      = "4.0",
      description  = "Python Trilinos Interface",
      author       = "Bill Spotz",
      author_email = "wfspotz@sandia.gov",
      package_dir  = {"PyTrilinos" : srcDir},
      packages     = ["PyTrilinos", "PyTrilinos.NOX"],
      ext_modules  = [ RawEpetra,    EpetraExt,     NOX_Epetra,
                       NOX_Abstract, NOX_Parameter, NOX_Solver,
                       NOX_StatusTest                          ]
      )
