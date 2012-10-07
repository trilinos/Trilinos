# -*- cmake -*-

# @HEADER
# ************************************************************************
#
#                PyTrilinos: Python Interface to Trilinos
#                   Copyright (2009) Sandia Corporation
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
# Questions? Contact Bill Spotz (wfspotz@sandia.gov)
#
# ************************************************************************
# @HEADER

# Include the package macros
INCLUDE(TribitsPackageMacros)

# Set the package name
TRIBITS_PACKAGE(PyTrilinos DISABLE_STRONG_WARNINGS)

IF(NOT BUILD_SHARED_LIBS)
  MESSAGE(FATAL_ERROR "PyTrilinos can only be built with shared libraries. Building of shared libraries is currently set to OFF. To enable shared libraries please set the cache variable \"BUILD_SHARED_LIBS\" to ON")
ENDIF()

# Set the package version number
SET(PyTrilinos_VERSION 4.8d)

# Define the list of PyTrilinos packages and modules.  Packages are
# the top-level names only, such as Epetra, AztecOO and NOX.  The list
# of modules is similar to the packages, except it contains the names
# of python modules that correspond to nested namespaces, such as
# NOX.__init__ or NOX.Epetra.Interface.
SET(PyTrilinos_PACKAGES "")
SET(PyTrilinos_MODULES  "")

IF(PyTrilinos_ENABLE_Teuchos)
  APPEND_SET(PyTrilinos_PACKAGES Teuchos)
  APPEND_SET(PyTrilinos_MODULES  Teuchos)
ENDIF(PyTrilinos_ENABLE_Teuchos)

IF(PyTrilinos_ENABLE_Epetra)
  APPEND_SET(PyTrilinos_PACKAGES Epetra)
  APPEND_SET(PyTrilinos_MODULES  Epetra)
ENDIF(PyTrilinos_ENABLE_Epetra)

IF(PyTrilinos_ENABLE_Triutils)
  APPEND_SET(PyTrilinos_PACKAGES TriUtils)
  APPEND_SET(PyTrilinos_MODULES  TriUtils)
ENDIF(PyTrilinos_ENABLE_Triutils)

IF(PyTrilinos_ENABLE_EpetraExt)
  APPEND_SET(PyTrilinos_PACKAGES EpetraExt)
  APPEND_SET(PyTrilinos_MODULES  EpetraExt)
ENDIF(PyTrilinos_ENABLE_EpetraExt)

IF(PyTrilinos_ENABLE_Isorropia)
  APPEND_SET(PyTrilinos_PACKAGES Isorropia)
  APPEND_SET(PyTrilinos_MODULES  Isorropia.__init__    )
  Append_SET(PyTrilinos_MODULES  Isorropia.NestedEpetra)
ENDIF(PyTrilinos_ENABLE_Isorropia)

IF(PyTrilinos_ENABLE_Pliris)
  APPEND_SET(PyTrilinos_PACKAGES Pliris)
  APPEND_SET(PyTrilinos_MODULES  Pliris)
ENDIF(PyTrilinos_ENABLE_Pliris)

IF(PyTrilinos_ENABLE_AztecOO)
  APPEND_SET(PyTrilinos_PACKAGES AztecOO)
  APPEND_SET(PyTrilinos_MODULES  AztecOO)
ENDIF(PyTrilinos_ENABLE_AztecOO)

IF(PyTrilinos_ENABLE_Galeri)
  APPEND_SET(PyTrilinos_PACKAGES Galeri)
  APPEND_SET(PyTrilinos_MODULES  Galeri)
ENDIF(PyTrilinos_ENABLE_Galeri)

IF(PyTrilinos_ENABLE_Amesos)
  APPEND_SET(PyTrilinos_PACKAGES Amesos)
  APPEND_SET(PyTrilinos_MODULES  Amesos)
ENDIF(PyTrilinos_ENABLE_Amesos)

IF(PyTrilinos_ENABLE_Ifpack)
  APPEND_SET(PyTrilinos_PACKAGES IFPACK)
  APPEND_SET(PyTrilinos_MODULES  IFPACK)
ENDIF(PyTrilinos_ENABLE_Ifpack)

IF(PyTrilinos_ENABLE_Komplex)
  APPEND_SET(PyTrilinos_PACKAGES Komplex)
  APPEND_SET(PyTrilinos_MODULES  Komplex)
ENDIF(PyTrilinos_ENABLE_Komplex)

IF(PyTrilinos_ENABLE_Anasazi)
  APPEND_SET(PyTrilinos_PACKAGES Anasazi)
  APPEND_SET(PyTrilinos_MODULES  Anasazi)
ENDIF(PyTrilinos_ENABLE_Anasazi)

IF(PyTrilinos_ENABLE_ML)
  APPEND_SET(PyTrilinos_PACKAGES ML)
  APPEND_SET(PyTrilinos_MODULES  ML)
ENDIF(PyTrilinos_ENABLE_ML)

IF(PyTrilinos_ENABLE_NOX)
  APPEND_SET(PyTrilinos_PACKAGES NOX)
  APPEND_SET(PyTrilinos_MODULES  NOX.__init__  )
  APPEND_SET(PyTrilinos_MODULES  NOX.Abstract  )
  APPEND_SET(PyTrilinos_MODULES  NOX.StatusTest)
  APPEND_SET(PyTrilinos_MODULES  NOX.Solver    )
  IF(NOX_ENABLE_Epetra)
    APPEND_SET(PyTrilinos_MODULES NOX.Epetra.__init__ )
    APPEND_SET(PyTrilinos_MODULES NOX.Epetra.Interface)
  ENDIF(NOX_ENABLE_Epetra)
ENDIF(PyTrilinos_ENABLE_NOX)

#IF(NOX_ENABLE_LOCA)
#  APPEND_SET(PyTrilinos_PACKAGES LOCA)
#  APPEND_SET(PyTrilinos_MODULES  LOCA.__init__         )
#  APPEND_SET(PyTrilinos_MODULES  LOCA.Abstract         )
#  APPEND_SET(PyTrilinos_MODULES  LOCA.Extended         )
#  APPEND_SET(PyTrilinos_MODULES  LOCA.BorderedSystem   )
#  APPEND_SET(PyTrilinos_MODULES  LOCA.Continuation     )
#  APPEND_SET(PyTrilinos_MODULES  LOCA.MultiContinuation)
#  APPEND_SET(PyTrilinos_MODULES  LOCA.Hopf             )
#  APPEND_SET(PyTrilinos_MODULES  LOCA.TimeDependent    )
#  APPEND_SET(PyTrilinos_MODULES  LOCA.Pitchfork        )
#  APPEND_SET(PyTrilinos_MODULES  LOCA.Homotopy         )
#  APPEND_SET(PyTrilinos_MODULES  LOCA.TurningPoint     )
#  #APPEND_SET(PyTrilinos_MODULES  LOCA.Bifurcation      )
#  #APPEND_SET(PyTrilinos_MODULES  LOCA.Chan             )
#  IF(NOX_ENABLE_Epetra)
#    APPEND_SET(PyTrilinos_MODULES  LOCA.Epetra.__init__ )
#    APPEND_SET(PyTrilinos_MODULES  LOCA.Epetra.Interface)
#  ENDIF(NOX_ENABLE_Epetra)
#ENDIF(NOX_ENABLE_LOCA)

# Python interpreter, library, and header files are required.
FIND_PACKAGE(PythonInterp 2.2 REQUIRED)
FIND_PACKAGE(PythonLibs 2.2 REQUIRED)
FIND_PACKAGE(PythonInclude REQUIRED)

# Python module numpy is required
FIND_PACKAGE(NumPy 1.0 REQUIRED)

# Simple wrapper interface generator (SWIG) is required
FIND_PACKAGE(SWIG 2.0.0 REQUIRED)
# ${SWIG_USE_FILE} is hard-coded to be the installed UseSWIG.cmake
# file.  We want to include the local, modified version.
INCLUDE(UseSWIG)

# Print the enabled PyTrilinos packages
MESSAGE(STATUS "Enabled PyTrilinos modules:")
MESSAGE(STATUS "    ${PyTrilinos_PACKAGES}")

# Determine the PyTrilinos default installation directory prefix.  If
# CMAKE_INSTALL_PREFIX is set by the user, then use it.  Else the
# prefix is set to the value of python's sys.prefix.
EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print sys.prefix"
  OUTPUT_VARIABLE PYTHON_PREFIX
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(PyTrilinos_DEFAULT_INSTALL_PREFIX ${PYTHON_PREFIX})
ELSE(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(PyTrilinos_DEFAULT_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# Set the PyTrilinos install prefix
SET(PyTrilinos_INSTALL_PREFIX ${PyTrilinos_DEFAULT_INSTALL_PREFIX}
  CACHE PATH "The path prefix for where PyTrilinos will be installed, e.g. /usr/local")

# Add the src directory, which builds the pytrilinos library,
# generates the wrapper code, builds the PyTrilinos python package and
# compiles the PyTrilinos extension modules
ADD_SUBDIRECTORY(src)

# Add the test scripts
TRIBITS_ADD_TEST_DIRECTORIES(test)

# Add the example scripts
TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_EXCLUDE_FILES(
  doc/PyTrilinos-SciProg/CoupledSolver.png
  doc/PyTrilinos-SciProg/CoupledSolver.pptx
  doc/PyTrilinos-SciProg/Response.pdf
  doc/PyTrilinos-SciProg/Response.tex
  doc/PyTrilinos-SciProg/cage12rcb8.png
  )

# Execute the package postprocessing
TRIBITS_PACKAGE_POSTPROCESS()
