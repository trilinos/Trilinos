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
INCLUDE(PackageMacros)

# Set the package name
PACKAGE(PyTrilinos DISABLE_STRONG_WARNINGS)

# Set the package version number
SET(${PACKAGE_NAME}_VERSION 4.2d)

# Define the list of PyTrilinos packages and modules.  Packages are
# the top-level names only, such as Epetra, AztecOO and NOX.  The list
# of modules is similar to the packages, except it contains the names
# of python modules that correspond to nested namespaces, such as
# NOX.__init__ or NOX.Epetra.
SET(${PACKAGE_NAME}_PACKAGES "")
SET(${PACKAGE_NAME}_MODULES  "")

IF(PyTrilinos_ENABLE_Teuchos)
  APPEND_SET(${PACKAGE_NAME}_PACKAGES Teuchos)
  APPEND_SET(${PACKAGE_NAME}_MODULES  Teuchos)
ENDIF(PyTrilinos_ENABLE_Teuchos)

IF(PyTrilinos_ENABLE_Epetra)
  APPEND_SET(${PACKAGE_NAME}_PACKAGES Epetra)
  APPEND_SET(${PACKAGE_NAME}_MODULES  Epetra)
ENDIF(PyTrilinos_ENABLE_Epetra)

IF(PyTrilinos_ENABLE_Triutils)
  APPEND_SET(${PACKAGE_NAME}_PACKAGES TriUtils)
  APPEND_SET(${PACKAGE_NAME}_MODULES  TriUtils)
ENDIF(PyTrilinos_ENABLE_Triutils)

IF(PyTrilinos_ENABLE_EpetraExt)
  APPEND_SET(${PACKAGE_NAME}_PACKAGES EpetraExt)
  APPEND_SET(${PACKAGE_NAME}_MODULES  EpetraExt)
ENDIF(PyTrilinos_ENABLE_EpetraExt)

IF(PyTrilinos_ENABLE_AztecOO)
  APPEND_SET(${PACKAGE_NAME}_PACKAGES AztecOO)
  APPEND_SET(${PACKAGE_NAME}_MODULES  AztecOO)
ENDIF(PyTrilinos_ENABLE_AztecOO)

IF(PyTrilinos_ENABLE_Galeri)
  APPEND_SET(${PACKAGE_NAME}_PACKAGES Galeri)
  APPEND_SET(${PACKAGE_NAME}_MODULES  Galeri)
ENDIF(PyTrilinos_ENABLE_Galeri)

IF(PyTrilinos_ENABLE_Amesos)
  APPEND_SET(${PACKAGE_NAME}_PACKAGES Amesos)
  APPEND_SET(${PACKAGE_NAME}_MODULES  Amesos)
ENDIF(PyTrilinos_ENABLE_Amesos)

IF(PyTrilinos_ENABLE_Ifpack)
  APPEND_SET(${PACKAGE_NAME}_PACKAGES IFPACK)
  APPEND_SET(${PACKAGE_NAME}_MODULES  IFPACK)
ENDIF(PyTrilinos_ENABLE_Ifpack)

IF(PyTrilinos_ENABLE_Komplex)
  APPEND_SET(${PACKAGE_NAME}_PACKAGES Komplex)
  APPEND_SET(${PACKAGE_NAME}_MODULES  Komplex)
ENDIF(PyTrilinos_ENABLE_Komplex)

IF(PyTrilinos_ENABLE_Anasazi)
  APPEND_SET(${PACKAGE_NAME}_PACKAGES Anasazi)
  APPEND_SET(${PACKAGE_NAME}_MODULES  Anasazi)
ENDIF(PyTrilinos_ENABLE_Anasazi)

IF(PyTrilinos_ENABLE_Pliris)
  APPEND_SET(${PACKAGE_NAME}_PACKAGES Pliris)
  APPEND_SET(${PACKAGE_NAME}_MODULES  Pliris)
ENDIF(PyTrilinos_ENABLE_Pliris)

IF(PyTrilinos_ENABLE_ML)
  APPEND_SET(${PACKAGE_NAME}_PACKAGES ML)
  APPEND_SET(${PACKAGE_NAME}_MODULES  ML)
ENDIF(PyTrilinos_ENABLE_ML)

IF(PyTrilinos_ENABLE_NOX)
  APPEND_SET(${PACKAGE_NAME}_PACKAGES NOX)
  APPEND_SET(${PACKAGE_NAME}_MODULES  NOX.__init__  )
  APPEND_SET(${PACKAGE_NAME}_MODULES  NOX.Abstract  )
  APPEND_SET(${PACKAGE_NAME}_MODULES  NOX.StatusTest)
  APPEND_SET(${PACKAGE_NAME}_MODULES  NOX.Solver    )
  IF(NOX_ENABLE_Epetra)
    APPEND_SET(${PACKAGE_NAME}_MODULES NOX.Epetra.__init__ )
    APPEND_SET(${PACKAGE_NAME}_MODULES NOX.Epetra.Interface)
  ENDIF(NOX_ENABLE_Epetra)
ENDIF(PyTrilinos_ENABLE_NOX)

IF(NOX_ENABLE_LOCA)
  APPEND_SET(${PACKAGE_NAME}_PACKAGES LOCA)
  APPEND_SET(${PACKAGE_NAME}_MODULES  LOCA.__init__         )
  #APPEND_SET(${PACKAGE_NAME}_MODULES  LOCA.Abstract         )
  #APPEND_SET(${PACKAGE_NAME}_MODULES  LOCA.Bifurcation      )
  #APPEND_SET(${PACKAGE_NAME}_MODULES  LOCA.BorderedSystem   )
  #APPEND_SET(${PACKAGE_NAME}_MODULES  LOCA.Chan             )
  #APPEND_SET(${PACKAGE_NAME}_MODULES  LOCA.Continuation     )
  #APPEND_SET(${PACKAGE_NAME}_MODULES  LOCA.Extended         )
  #APPEND_SET(${PACKAGE_NAME}_MODULES  LOCA.Homotopy         )
  #APPEND_SET(${PACKAGE_NAME}_MODULES  LOCA.Hopf             )
  #APPEND_SET(${PACKAGE_NAME}_MODULES  LOCA.MultiContinuation)
  #APPEND_SET(${PACKAGE_NAME}_MODULES  LOCA.Pitchfork        )
  #APPEND_SET(${PACKAGE_NAME}_MODULES  LOCA.TimeDependent    )
  #APPEND_SET(${PACKAGE_NAME}_MODULES  LOCA.TurningPoint     )
  #IF(NOX_ENABLE_Epetra)
  #  APPEND_SET(${PACKAGE_NAME}_MODULES  LOCA.Epetra.__init__ )
  #  APPEND_SET(${PACKAGE_NAME}_MODULES  LOCA.Epetra.Interface)
  #ENDIF(NOX_ENABLE_Epetra)
ENDIF(NOX_ENABLE_LOCA)

# Print the enabled PyTrilinos packages
MESSAGE(STATUS "Enabled ${PACKAGE_NAME} modules:")
MESSAGE(STATUS "    ${${PACKAGE_NAME}_PACKAGES}")

# Python interpreter and header files are required
FIND_PACKAGE(PythonInterp 2.2 REQUIRED)
FIND_PACKAGE(PythonInclude REQUIRED)

# Python module numpy is required
FIND_PACKAGE(NumPy 1.0 REQUIRED)

# Simple wrapper interface generator (SWIG) is required
FIND_PACKAGE(SWIG 1.3.37 REQUIRED)
# ${SWIG_USE_FILE} is hard-coded to be the installed UseSWIG.cmake
# file.  We want to include the local, modified version.
INCLUDE(UseSWIG)

# Add the src directory, which builds the pytrilinos library,
# generates the wrapper code, builds the PyTrilinos python package and
# compiles the PyTrilinos extension modules
ADD_SUBDIRECTORY(src)

# Add the test scripts
PACKAGE_ADD_TEST_DIRECTORIES(test example)

# Add the example scripts
#PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

# Execute the package postprocessing
PACKAGE_POSTPROCESS()
