# -*- cmake -*-

# @HEADER
# ***********************************************************************
#
#          PyTrilinos: Python Interfaces to Trilinos Packages
#                 Copyright (2014) Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia
# Corporation, the U.S. Government retains certain rights in this
# software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact William F. Spotz (wfspotz@sandia.gov)
#
# ***********************************************************************
# @HEADER

# Include the package macros

# Set the package name
TRIBITS_PACKAGE(PyTrilinos DISABLE_STRONG_WARNINGS)

IF(NOT BUILD_SHARED_LIBS)
  MESSAGE(FATAL_ERROR "PyTrilinos can only be built with shared libraries. Building of shared libraries is currently set to OFF. To enable shared libraries please set the cache variable \"BUILD_SHARED_LIBS\" to ON")
ENDIF()

# Set the package version number
SET(PyTrilinos_VERSION ${Trilinos_VERSION})

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

# Define the list of PyTrilinos packages and modules.  Packages are
# the top-level names only, such as Epetra, AztecOO and NOX.  The list
# of modules is similar to the packages, except it contains the names
# of python modules that correspond to nested namespaces, such as
# NOX.__init__ or NOX.Epetra.Interface.
SET(PyTrilinos_PACKAGES "")
SET(PyTrilinos_MODULES  "")

APPEND_SET(PyTrilinos_PACKAGES Teuchos)
APPEND_SET(PyTrilinos_MODULES  Teuchos)

IF(PyTrilinos_ENABLE_Epetra)
  APPEND_SET(PyTrilinos_PACKAGES Epetra)
  APPEND_SET(PyTrilinos_MODULES  Epetra)
ENDIF()

IF(PyTrilinos_ENABLE_Triutils)
  APPEND_SET(PyTrilinos_PACKAGES TriUtils)
  APPEND_SET(PyTrilinos_MODULES  TriUtils)
ENDIF()

IF(PyTrilinos_ENABLE_Tpetra)
  APPEND_SET(PyTrilinos_PACKAGES Tpetra)
  APPEND_SET(PyTrilinos_MODULES  Tpetra)
ENDIF()

IF(PyTrilinos_ENABLE_Domi)
  APPEND_SET(PyTrilinos_PACKAGES Domi)
  APPEND_SET(PyTrilinos_MODULES  Domi)
ENDIF()

IF(PyTrilinos_ENABLE_EpetraExt)
  APPEND_SET(PyTrilinos_PACKAGES EpetraExt)
  APPEND_SET(PyTrilinos_MODULES  EpetraExt)
ENDIF()

IF(PyTrilinos_ENABLE_Isorropia)
  APPEND_SET(PyTrilinos_PACKAGES Isorropia)
  APPEND_SET(PyTrilinos_MODULES  Isorropia)
  IF(PyTrilinos_ENABLE_Epetra)
    Append_SET(PyTrilinos_MODULES Isorropia.Epetra)
  ENDIF()
ENDIF()

IF(PyTrilinos_ENABLE_Pliris)
  APPEND_SET(PyTrilinos_PACKAGES Pliris)
  APPEND_SET(PyTrilinos_MODULES  Pliris)
ENDIF()

IF(PyTrilinos_ENABLE_AztecOO)
  APPEND_SET(PyTrilinos_PACKAGES AztecOO)
  APPEND_SET(PyTrilinos_MODULES  AztecOO)
ENDIF()

IF(PyTrilinos_ENABLE_Galeri)
  APPEND_SET(PyTrilinos_PACKAGES Galeri)
  APPEND_SET(PyTrilinos_MODULES  Galeri)
ENDIF()

IF(PyTrilinos_ENABLE_Amesos)
  APPEND_SET(PyTrilinos_PACKAGES Amesos)
  APPEND_SET(PyTrilinos_MODULES  Amesos)
ENDIF()

IF(PyTrilinos_ENABLE_Ifpack)
  APPEND_SET(PyTrilinos_PACKAGES IFPACK)
  APPEND_SET(PyTrilinos_MODULES  IFPACK)
ENDIF()

IF(PyTrilinos_ENABLE_Komplex)
  APPEND_SET(PyTrilinos_PACKAGES Komplex)
  APPEND_SET(PyTrilinos_MODULES  Komplex)
ENDIF()

IF(PyTrilinos_ENABLE_Anasazi)
  APPEND_SET(PyTrilinos_PACKAGES Anasazi)
  APPEND_SET(PyTrilinos_MODULES  Anasazi)
ENDIF()

IF(PyTrilinos_ENABLE_ML)
  APPEND_SET(PyTrilinos_PACKAGES ML)
  APPEND_SET(PyTrilinos_MODULES  ML)
ENDIF()

IF(PyTrilinos_ENABLE_NOX)
  APPEND_SET(PyTrilinos_PACKAGES NOX)
  APPEND_SET(PyTrilinos_MODULES  NOX           )
  APPEND_SET(PyTrilinos_MODULES  NOX.Abstract  )
  APPEND_SET(PyTrilinos_MODULES  NOX.StatusTest)
  APPEND_SET(PyTrilinos_MODULES  NOX.Solver    )
  ASSERT_DEFINED(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA)
  IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA)
    APPEND_SET(PyTrilinos_MODULES NOX.Epetra          )
    APPEND_SET(PyTrilinos_MODULES NOX.Epetra.Interface)
  ENDIF()
  IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_PETSC)
    APPEND_SET(PyTrilinos_MODULES NOX.PETSc          )
    APPEND_SET(PyTrilinos_MODULES NOX.PETSc.Interface)
  ENDIF()
ENDIF()

IF(NOX_ENABLE_LOCA)
  APPEND_SET(PyTrilinos_PACKAGES LOCA                                )
  APPEND_SET(PyTrilinos_MODULES  LOCA                                )
  APPEND_SET(PyTrilinos_MODULES  LOCA.Abstract                       )
  APPEND_SET(PyTrilinos_MODULES  LOCA.AnasaziOperator                )
  APPEND_SET(PyTrilinos_MODULES  LOCA.Bifurcation                    )
  APPEND_SET(PyTrilinos_MODULES  LOCA.BorderedSolver                 )
  APPEND_SET(PyTrilinos_MODULES  LOCA.BorderedSystem                 )
  APPEND_SET(PyTrilinos_MODULES  LOCA.Eigensolver                    )
  APPEND_SET(PyTrilinos_MODULES  LOCA.EigenvalueSort                 )
  APPEND_SET(PyTrilinos_MODULES  LOCA.Extended                       )
  APPEND_SET(PyTrilinos_MODULES  LOCA.Homotopy                       )
  APPEND_SET(PyTrilinos_MODULES  LOCA.Hopf                           )
  APPEND_SET(PyTrilinos_MODULES  LOCA.Hopf.MinimallyAugmented        )
  APPEND_SET(PyTrilinos_MODULES  LOCA.Hopf.MooreSpence               )
  APPEND_SET(PyTrilinos_MODULES  LOCA.MultiContinuation              )
  APPEND_SET(PyTrilinos_MODULES  LOCA.MultiPredictor                 )
  APPEND_SET(PyTrilinos_MODULES  LOCA.Parameter                      )
  APPEND_SET(PyTrilinos_MODULES  LOCA.PhaseTransition                )
  APPEND_SET(PyTrilinos_MODULES  LOCA.Pitchfork                      )
  APPEND_SET(PyTrilinos_MODULES  LOCA.Pitchfork.MinimallyAugmented   )
  APPEND_SET(PyTrilinos_MODULES  LOCA.Pitchfork.MooreSpence          )
  APPEND_SET(PyTrilinos_MODULES  LOCA.SaveEigenData                  )
  APPEND_SET(PyTrilinos_MODULES  LOCA.StatusTest                     )
  APPEND_SET(PyTrilinos_MODULES  LOCA.StepSize                       )
  APPEND_SET(PyTrilinos_MODULES  LOCA.TimeDependent                  )
  APPEND_SET(PyTrilinos_MODULES  LOCA.TurningPoint                   )
  APPEND_SET(PyTrilinos_MODULES  LOCA.TurningPoint.MinimallyAugmented)
  APPEND_SET(PyTrilinos_MODULES  LOCA.TurningPoint.MooreSpence       )
  IF(NOX_ENABLE_Epetra)
    APPEND_SET(PyTrilinos_MODULES  LOCA.Epetra          )
    APPEND_SET(PyTrilinos_MODULES  LOCA.Epetra.Interface)
  ENDIF()
ENDIF()

# Python interpreter and header files are required
FIND_PACKAGE(PythonInterp 2.5 REQUIRED)
FIND_PACKAGE(PythonInclude REQUIRED)

# Python module numpy is required
FIND_PACKAGE(NumPy 1.0 REQUIRED)

# Simple wrapper interface generator (SWIG) is required. SWIG 3.0.11
# or higher is required to support the %module(moduleimport=...)
# option
SET(SWIG_VERSION_FLOOR 3.0.11)
FIND_PACKAGE(SWIG REQUIRED)
IF(${SWIG_VERSION} VERSION_LESS ${SWIG_VERSION_FLOOR})
  MESSAGE(FATAL_ERROR
    "SWIG version " ${SWIG_VERSION}
    " is less than minimum version " ${SWIG_VERSION_FLOOR}
    )
ENDIF()
MESSAGE(STATUS "Found swig: " ${SWIG_VERSION})
# ${SWIG_USE_FILE} is hard-coded to be the installed UseSWIG.cmake
# file.  We want to include the local, modified version.
INCLUDE(UseSWIG)

# Python module mpi4py is optional, and only of interest if MPI is
# enabled
IF(TPL_ENABLE_MPI)
  SET(TPL_ENABLE_MPI4PY ON CACHE BOOL "Use Mpi4Py")
  IF(TPL_ENABLE_MPI4PY)
    FIND_PACKAGE(Mpi4Py)
  ELSE()
    SET(Mpi4Py_FOUND OFF)
  ENDIF()
ENDIF()

# Python module petsc4py is required if NOX_ENABLE_PETSC is true
IF(NOX_ENABLE_PETSC)
  FIND_PACKAGE(PETSc4Py REQUIRED)
ENDIF()

# Print the enabled PyTrilinos packages
MESSAGE(STATUS "Enabled PyTrilinos packages:")
MESSAGE(STATUS "    ${PyTrilinos_PACKAGES}")

# Determine the PyTrilinos default installation directory prefix.  If
# CMAKE_INSTALL_PREFIX is set by the user, then use it.  Else the
# prefix is set to the value of python's sys.prefix.
EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print(sys.prefix)"
  OUTPUT_VARIABLE PYTHON_PREFIX
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(PyTrilinos_DEFAULT_INSTALL_PREFIX ${PYTHON_PREFIX})
ELSE()
  SET(PyTrilinos_DEFAULT_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
ENDIF()

# Set the PyTrilinos install prefix
SET(PyTrilinos_INSTALL_PREFIX ${PyTrilinos_DEFAULT_INSTALL_PREFIX}
  CACHE PATH "The path prefix for where PyTrilinos will be installed, e.g. /usr/local")

# The default status for whether we build docstrings for PyTrilinos is
# based on whether the user has doxygen installed.
FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
  IF(${DOXYGEN_VERSION} VERSION_EQUAL 1.8.13 AND PyTrilinos_DOCSTRINGS)
    MESSAGE(FATAL_ERROR 
      "Doxygen version ${DOXYGEN_VERSION} has a bug that makes it incompatible with "
      "PyTrilinos. To work around this, you may disable PyTrilinos docstrings by setting "
      "PyTrilinos_DOCSTRINGS=OFF. "
      "Alternatively, you may choose to use a different version of doxygen. If you "
      "choose this route, be sure that the new version of doxygen is found first in "
      "your environment path.")
  ELSE()
    SET(PyTrilinos_DOCSTRINGS_DEFAULT ON)
  ENDIF()
ELSE()
  SET(PyTrilinos_DOCSTRINGS_DEFAULT OFF)
  MESSAGE(STATUS "doxygen not found!")
ENDIF()

# Set the PyTrilinos docstrings flag
SET(PyTrilinos_DOCSTRINGS ${PyTrilinos_DOCSTRINGS_DEFAULT} CACHE BOOL
    "a flag that determines whether Python docstrings will be built using doxygen")
IF(PyTrilinos_DOCSTRINGS AND NOT DOXYGEN_FOUND)
  MESSAGE(ERROR "PyTrilinos docstrings enabled but doxygen not found")
ENDIF()

# Add the util directory
ADD_SUBDIRECTORY(util)

# Add the doc directory, which builds the python documentation strings
IF(SWIG_VERSION VERSION_LESS 4.0)
  ADD_SUBDIRECTORY(doc)
ENDIF()

# Add the src directory, which builds the pytrilinos library,
# generates the wrapper code, builds the PyTrilinos python package and
# compiles the PyTrilinos extension modules
ADD_SUBDIRECTORY(src)

# Add the test scripts
TRIBITS_ADD_TEST_DIRECTORIES(test)

# Add the example scripts
TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

# Execute the package postprocessing
TRIBITS_PACKAGE_POSTPROCESS()
