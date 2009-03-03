
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
SET(DEFINES    "")
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
SET(EXTRA_LINK_ARGS "")
IF(APPLE)
  IF(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    APPEND_SET(EXTRA_LINK_ARGS "-undefined dynamic_lookup")
  ENDIF()
ENDIF(APPLE)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

######################################
# Core PyTrilinos headers and sources
######################################

APPEND_SET(HEADERS
  PythonException.h
  FILEstream.h
  NumPyImporter.h
  numpy_include.h
  )

APPEND_SET(SOURCES
  PythonException.cpp
  FILEstream.cpp
  )

#########################################
# PyTrilinos+Teuchos headers and sources
#########################################

IF(Trilinos_ENABLE_Teuchos)

  APPEND_SET(HEADERS
    Teuchos_PythonParameter.h
    )

  APPEND_SET(SOURCES
    Teuchos_PythonParameter.cpp
    )

  APPEND_SET(DEFINES
    -DHAVE_TEUCHOS
    )

ENDIF(Trilinos_ENABLE_Teuchos)

########################################
# PyTrilinos+Epetra headers and sources
########################################

IF(Trilinos_ENABLE_Epetra)

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

  APPEND_SET(DEFINES
    -DHAVE_EPETRA
    )

ENDIF(Trilinos_ENABLE_Epetra)

####################################
# PyTrilinos+ML headers and sources
####################################

IF(Trilinos_ENABLE_ML)

  APPEND_SET(HEADERS
    MLAPI_PyMatrix.h
    )

ENDIF(Trilinos_ENABLE_ML)

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
  PROPERTIES LINK_FLAGS "${EXTRA_LINK_ARGS}"
  )

#
# Define the PyTrilinos swig setup
SET(CMAKE_SWIG_FLAGS "")
GET_FILENAME_COMPONENT(PYTRILINOS_DOXDIR
  ${CMAKE_CURRENT_SOURCE_DIR}/../doc/Doxygen ABSOLUTE
  )
APPEND_SET(CMAKE_SWIG_FLAGS -I${PYTRILINOS_DOXDIR})
SET(CMAKE_SWIG_OUTDIR PyTrilinos)

#
# Define the swig-generated PyTrilinos.Teuchos module
IF(Trilinos_ENABLE_Teuchos)

  SET_SOURCE_FILES_PROPERTIES(Teuchos.i
    PROPERTIES CPLUSPLUS ON
    SWIG_FLAGS "${DEFINES}"
    )
  SWIG_ADD_MODULE(Teuchos python Teuchos.i)
  SWIG_LINK_LIBRARIES(Teuchos ${TEUCHOS_LIBRARIES})
  SWIG_LINK_LIBRARIES(Teuchos pytrilinos)
  SET_TARGET_PROPERTIES(${SWIG_MODULE_Teuchos_REAL_NAME}
    PROPERTIES LIBRARY_OUTPUT_DIRECTORY PyTrilinos
    )
  SET_TARGET_PROPERTIES(${SWIG_MODULE_Teuchos_REAL_NAME}
    PROPERTIES LINK_FLAGS "${EXTRA_LINK_ARGS}"
    )

ENDIF(Trilinos_ENABLE_Teuchos)

#
# Define the swig-generated PyTrilinos.Epetra module
IF(Trilinos_ENABLE_Epetra)

  SET_SOURCE_FILES_PROPERTIES(Epetra.i
    PROPERTIES CPLUSPLUS ON
    SWIG_FLAGS "-DHAVE_EPETRA"
    )
  SWIG_ADD_MODULE(Epetra python Epetra.i)
  SWIG_LINK_LIBRARIES(Epetra ${EPETRA_LIBRARIES})
  SWIG_LINK_LIBRARIES(Epetra pytrilinos)
  SET_TARGET_PROPERTIES(${SWIG_MODULE_Epetra_REAL_NAME}
    PROPERTIES LIBRARY_OUTPUT_DIRECTORY PyTrilinos
    )
  SET_TARGET_PROPERTIES(${SWIG_MODULE_Epetra_REAL_NAME}
    PROPERTIES LINK_FLAGS "${EXTRA_LINK_ARGS}"
    )

ENDIF(Trilinos_ENABLE_Epetra)

##############

ADD_SUBDIRECTORY(PyTrilinos)
