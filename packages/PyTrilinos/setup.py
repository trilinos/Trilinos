#! /usr/bin/env python

# System imports
from   distutils.core import *
from   distutils      import sysconfig
import commands
import os
import string
import sys

# Build the python library directory name and library name.  These will be used
# when the linker is invoked to create the python extensions.
pythonDir = [sysconfig.get_config_var('LIBPL'  )      ]
pythonLib = [sysconfig.get_config_var('LIBRARY')[3:-2]]

# Any information that needs to be transferred from the autotooled Makefile is
# written to file setup.txt using python syntax to define a dictionary.  The
# keys of this 'makeInfo' dictionary are variable names and the corresponding
# values represent the data that will be needed by this setup.py script.
f = open("setup.txt")
makeInfo = f.readlines()
f.close()
makeInfo = eval(string.join(makeInfo))

# Certain directory paths are needed by setup.py.  pyTDir is the path for the
# PyTrilinos package, pakDir is the path for the Trilinos package directory,
# srcDir is the PyTrilinos source directory path, and noxDir is the PyTrilinos
# NOX source directory path.
pyTDir = makeInfo['srcdir']
pakDir = os.path.split(pyTDir      )[0]
srcDir = os.path.join( pyTDir,"src")
noxDir = os.path.join( srcDir,"NOX")

# Define the include paths required by various packages.  Each of these is
# defined as a list of a single or more strings.  Thus they can be added
# together to yield a list of multiple strings.
epetraInc    = [os.path.join(pakDir, "epetra",    "src"             )  ]
epetraExtInc = [os.path.join(pakDir, "epetraext", "src", "transform"), \
                os.path.join(pakDir, "epetraext", "src", "coloring" )  ]
noxInc       = [os.path.join(pakDir, "nox",       "src"             )  ]
noxEpetraInc = [os.path.join(pakDir, "nox",       "src-epetra"      )  ]

# Define the library search directories needed to link to various package
# libraries.  Each of these is defined as a list of a single string.  Thus they
# can be added together to yield a list of multiple strings.
epetraLibDir    = [os.path.join("..", "epetra",    "src"       )]
epetraExtLibDir = [os.path.join("..", "epetraext", "src"       )]
noxLibDir       = [os.path.join("..", "nox",       "src"       )]
noxEpetraLibDir = [os.path.join("..", "nox",       "src-epetra")]
aztecLibDir     = [os.path.join("..", "aztecoo",   "src"       )]
ifpackLibDir    = [os.path.join("..", "ifpack",    "src"       )]

# Define the library names for various packages.  Each of these is defined as a
# list of a single string.  Thus they can be added together to yield a list of
# multiple strings.
epetraLib    = ["epetra"   ]
epetraExtLib = ["epetraext"]
noxEpetraLib = ["noxepetra"]
noxLib       = ["nox"      ]
aztecLib     = ["aztecoo"  ]
ifpackLib    = ["ifpack"   ]

# Obtain the directory path and library name for the SWIG runtime library.
# Currently, this requires SWIG Version 1.3.21, using flags that are depricated
# in SWIG Version 1.3.22.  Ultimately, this will be replaced with a locally
# generated runtime library, but that will require libtool, which is currently
# beyond the scope of the project.
swigStrings = commands.getoutput("swig -python -ldflags").split()
swigDir     = [swigStrings[0][2:]]
swigLib     = [swigStrings[1][2:]]
sysName = os.uname()[0]
if sysName == "Linux":
    swigPath = ["-Wl,-rpath"] + swigDir
else:
    swigPath = [ ]

# Standard libraries.  This is currently a hack.  The library "stdc++" is added
# to the standard library list for a case where we know it needs it.  Also, the
# SWIG runtime library is added to the list.
if sysName == "Linux":
    stdCXX = ["stdc++"]
else:
    stdCXX = [ ]
stdLibs    = stdCXX + swigLib

# Create the extra arguments list and complete the standard liraries list.  This
# is accomplished by looping over the arguments in FLIB and adding them to the
# appropriate list.
extraArgs = swigPath
flibs = makeInfo['FLIBS'].split()
for lib in flibs:
    if lib[:2] == "-l":
        stdLibs.append(lib[2:])
    else:
        extraArgs.append(lib)

# Define the strings that refer to the specified source files.  These should be
# combined together as needed in a list as the second argument to the Extension
# constructor.
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
                      libraries       = stdLibs + epetraLib + pythonLib,
                      extra_link_args = extraArgs
                      )

# EpetraExt extension module
EpetraExt = Extension("PyTrilinos._EpetraExt",
                      [epetraExtWrap],
                      include_dirs    = epetraInc + epetraExtInc,
                      library_dirs    = epetraLibDir + epetraExtLibDir + swigDir + \
                                        pythonDir,
                      libraries       = stdLibs + epetraLib + epetraExtLib + \
                                        pythonLib,
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
                       include_dirs    = epetraInc + noxInc + noxEpetraInc,
                       library_dirs    = epetraLibDir + noxLibDir + noxEpetraLibDir + \
                                         swigDir,
                       libraries       = stdLibs + epetraLib + noxLib + noxEpetraLib + \
                                         aztecLib + ifpackLib,
                       extra_link_args = extraArgs
                       )

# NOX.Abstract extension module
NOX_Abstract = Extension("PyTrilinos.NOX._Abstract",
                         [noxAbstractWrap],
                         include_dirs    = noxInc,
                         library_dirs    = noxLibDir + swigDir,
                         libraries       = stdLibs + noxLib,
                         extra_link_args = extraArgs
                         )

# NOX.Parameter extension module
NOX_Parameter = Extension("PyTrilinos.NOX._Parameter",
                          [noxParameterWrap],
                          include_dirs    = noxInc,
                          library_dirs    = noxLibDir + swigDir,
                          libraries       = stdLibs + noxLib,
                          extra_link_args = extraArgs
                          )

# NOX.Solver extension module
NOX_Solver = Extension("PyTrilinos.NOX._Solver",
                       [noxSolverWrap],
                       include_dirs    = noxInc,
                       library_dirs    = noxLibDir + swigDir,
                       libraries       = stdLibs + noxLib,
                       extra_link_args = extraArgs
                       )

# NOX.StatusTest extension module
NOX_StatusTest = Extension("PyTrilinos.NOX._StatusTest",
                           [noxStatusTestWrap],
                           include_dirs    = noxInc,
                           library_dirs    = noxLibDir + swigDir,
                           libraries       = stdLibs + noxLib,
                           extra_link_args = extraArgs
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
