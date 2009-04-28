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

INCLUDE(PackageLibraryMacros)

#
# Package-specific configuration options
IF(Trilinos_ENABLE_Epetra)
  SET(HAVE_EPETRA ON)
ENDIF(Trilinos_ENABLE_Epetra)

IF(Trilinos_ENABLE_Teuchos)
  SET(HAVE_TEUCHOS ON)
ENDIF(Trilinos_ENABLE_Teuchos)

IF(NOX_ENABLE_Epetra)
  SET(HAVE_NOX_EPETRA ON)
ENDIF(NOX_ENABLE_Epetra)

IF(TPL_ENABLE_MPI)
  SET(HAVE_MPI ON)
ENDIF(TPL_ENABLE_MPI)

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# On Mac OS X Gnu compilers, add dynamic lookup for undefined symbols
# to the pytrilinos library and PyTrilinos extension modules
SET(EXTRA_LINK_ARGS "")
IF(APPLE)
  IF(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    APPEND_SET(EXTRA_LINK_ARGS "-undefined dynamic_lookup")
  ENDIF()
ENDIF(APPLE)

#
# Add the current binary and source directories, the python include
# path, and the NumPy include path to the include directories
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${NumPy_INCLUDE_DIR})

#
# Copy the pyfragments.swg file to the build directory
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/pyfragments.swg.in
  ${CMAKE_CURRENT_BINARY_DIR}/pyfragments.swg)

#
# Initialize the headers and sources
SET(HEADERS "")
SET(SOURCES "")
SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

#
# Core PyTrilinos headers and sources
APPEND_SET(HEADERS
  NumPyImporter.h
  PythonException.h
  FILEstream.h
  numpy_include.h
  )

APPEND_SET(SOURCES
  NumPyImporter.cpp
  PythonException.cpp
  FILEstream.cpp
  )

#
# PyTrilinos + Teuchos headers and sources
IF(PyTrilinos_ENABLE_Teuchos)

  APPEND_SET(HEADERS
    Teuchos_PythonParameter.h
    )

  APPEND_SET(SOURCES
    Teuchos_PythonParameter.cpp
    )

ENDIF(PyTrilinos_ENABLE_Teuchos)

#
# PyTrilinos + Epetra headers and sources
IF(PyTrilinos_ENABLE_Epetra)

  APPEND_SET(HEADERS
    Epetra_NumPyFEVector.h
    Epetra_NumPyIntSerialDenseMatrix.h
    Epetra_NumPyIntSerialDenseVector.h
    Epetra_NumPyIntVector.h
    Epetra_NumPyMultiVector.h
    Epetra_NumPySerialDenseMatrix.h
    Epetra_NumPySerialDenseVector.h
    Epetra_NumPyVector.h
    )

  APPEND_SET(SOURCES
    Epetra_NumPyFEVector.cpp
    Epetra_NumPyIntSerialDenseMatrix.cpp
    Epetra_NumPyIntSerialDenseVector.cpp
    Epetra_NumPyIntVector.cpp
    Epetra_NumPyMultiVector.cpp
    Epetra_NumPySerialDenseMatrix.cpp
    Epetra_NumPySerialDenseVector.cpp
    Epetra_NumPyVector.cpp
    )

ENDIF(PyTrilinos_ENABLE_Epetra)

#
# PyTrilinos + ML headers and sources
IF(PyTrilinos_ENABLE_ML)

  APPEND_SET(HEADERS
    MLAPI_PyMatrix.h
    )

ENDIF(PyTrilinos_ENABLE_ML)

#
# Define the target for the swig runtime header
SET(SWIG_RUNTIME swigpyrun.h)
ADD_CUSTOM_COMMAND(
  OUTPUT ${SWIG_RUNTIME}
  COMMAND ${SWIG_EXECUTABLE} -python -external-runtime
  )

#
# Define the targets for the PyTrilinos library
PACKAGE_ADD_LIBRARY(
  pytrilinos
  HEADERS ${HEADERS} ${SWIG_RUNTIME}
  SOURCES ${SOURCES}
  )

SET_TARGET_PROPERTIES(pytrilinos
  PROPERTIES LINK_FLAGS "${EXTRA_LINK_ARGS}"
  )

#
# Define the PyTrilinos swig setup
SET(CMAKE_SWIG_FLAGS "")
GET_FILENAME_COMPONENT(PYTRILINOS_DOXDIR
  ${CMAKE_CURRENT_SOURCE_DIR}/../doc/Doxygen ABSOLUTE
  )
APPEND_SET(CMAKE_SWIG_FLAGS -I${PYTRILINOS_DOXDIR})

#
# Loop over the PyTrilinos-enabled Trilinos modules and define the
# swig-generated PyTrilinos extension modules
FOREACH(MODULE ${PyTrilinos_MODULES})
  # Alternate versions of the module name
  STRING(TOUPPER ${MODULE} MODULE_UPPER)
  STRING(REPLACE "." "_" MODULE_NODOT ${MODULE})
  SET_SOURCE_FILES_PROPERTIES(${MODULE}.i PROPERTIES
    CPLUSPLUS ON
    )
  # Determine the swig module output directory, name and unique target name
  SWIG_MODULE_GET_OUTDIR_AND_MODULE(${MODULE}.i MODULE_OUTDIR MODULE_NAME)
  IF(TARGET ${MODULE_NAME})
    SET(PyTrilinos_TARGET_NAME ${MODULE_NODOT})
  ELSE(TARGET ${MODULE_NAME})
    SET(PyTrilinos_TARGET_NAME ${MODULE_NAME})
  ENDIF(TARGET ${MODULE_NAME})
  # Add the SWIG module
  SWIG_ADD_MODULE(${PyTrilinos_TARGET_NAME} python "${MODULE_OUTDIR}"
    "${MODULE_NAME}" ${MODULE}.i)
  SWIG_LINK_LIBRARIES(${PyTrilinos_TARGET_NAME} pytrilinos)
  # Choosing the Trilinos libraries to link against requires a little bit
  # of logic
  IF("${MODULE_OUTDIR}" MATCHES "PyTrilinos/NOX/Epetra")
    SET(TRILINOS_LIBS ${NOXEPETRA_LIBRARIES})
  ELSEIF("${MODULE_OUTDIR}" MATCHES "PyTrilinos/NOX")
    SET(TRILINOS_LIBS ${NOX_LIBRARIES})
  ELSEIF("${MODULE_OUTDIR}" MATCHES "PyTrilinos/LOCA")
    SET(TRILINOS_LIBS ${LOCA_LIBRARIES})
  ELSE("${MODULE_OUTDIR}")
    SET(TRILINOS_LIBS ${${MODULE_UPPER}_LIBRARIES})
  ENDIF("${MODULE_OUTDIR}" MATCHES "PyTrilinos/NOX/Epetra")
  SWIG_LINK_LIBRARIES(${PyTrilinos_TARGET_NAME} ${TRILINOS_LIBS})
  SET_TARGET_PROPERTIES(${PyTrilinos_TARGET_NAME} PROPERTIES
    LINK_FLAGS "${EXTRA_LINK_ARGS}"
    )
ENDFOREACH(MODULE ${PyTrilinos_MODULES})

#
# Add the PyTrilinos subdirectory, which is top-level package
# directory
ADD_SUBDIRECTORY(PyTrilinos)
