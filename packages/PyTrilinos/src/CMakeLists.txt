
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS    "")
SET(SOURCES    "")
SET(SWIG_OPTS  "")
SET(INTERFACES "")

# Add the current binary directory to the include directories
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# Add the python include path to the include directories
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

# Add the NumPy include path to the include directories
INCLUDE_DIRECTORIES(${NumPy_INCLUDE_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

# On Mac OS X Gnu compilers, add dynamic lookup for undefined symbols
# to the pytrilinos library
SET(UNDEFINED_DYNAMIC_LOOKUP "")
IF(APPLE)
  IF(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    APPEND_SET(UNDEFINED_DYNAMIC_LOOKUP "-undefined dynamic_lookup")
  ENDIF()
ENDIF(APPLE)


#
# Core headers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Epetra_NumPyFEVector.h
  Epetra_NumPyIntSerialDenseMatrix.h
  Epetra_NumPyIntSerialDenseVector.h
  Epetra_NumPyIntVector.h
  Epetra_NumPyMultiVector.h
  Epetra_NumPySerialDenseMatrix.h
  Epetra_NumPySerialDenseVector.h
  Epetra_NumPyVector.h
  FILEstream.h
  MLAPI_PyMatrix.h
  NumPyImporter.h
  PythonException.h
  Teuchos_PythonParameter.h
  numpy_include.h
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
  FILEstream.cpp
  PythonException.cpp
  Teuchos_PythonParameter.cpp
  )

######################################

#APPEND_SET(HEADERS
#  )

#APPEND_SET(SOURCES
#  )

######################################
IF (TPL_ENABLE_MPI)
ENDIF()

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
  PROPERTIES LINK_FLAGS "${PYTHON_LIBRARIES}"
  )

#
# Define the PyTrilinos swig setup
SET(CMAKE_SWIG_FLAGS "")
GET_FILENAME_COMPONENT(PYTRILINOS_DOXDIR
  ${CMAKE_CURRENT_SOURCE_DIR}/../doc/Doxygen ABSOLUTE
)
APPEND_SET(CMAKE_SWIG_FLAGS -I${PYTRILINOS_DOXDIR})

#
# Define the swig-generated PyTrilinos.Epetra module
SET_SOURCE_FILES_PROPERTIES(Epetra.i
  PROPERTIES CPLUSPLUS ON
             SWIG_FLAGS "${CMAKE_SWIG_FLAGS}"
  )
SWIG_ADD_MODULE(Epetra python Epetra.i)
SWIG_LINK_LIBRARIES(Epetra ${EPETRA_LIBRARIES})
SWIG_LINK_LIBRARIES(Epetra pytrilinos)
SWIG_LINK_LIBRARIES(Epetra ${PYTHON_LIBRARIES})
