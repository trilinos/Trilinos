
INCLUDE(PackageLibraryMacros)

#
# Package-specific configuration options
IF(Trilinos_ENABLE_Epetra)
  SET(HAVE_EPETRA ON)
ENDIF(Trilinos_ENABLE_Epetra)

IF(Trilinos_ENABLE_Teuchos)
  SET(HAVE_TEUCHOS ON)
ENDIF(Trilinos_ENABLE_Teuchos)

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
# Initialize the headers and sources
SET(HEADERS "")
SET(SOURCES "")
SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

#
# Core PyTrilinos headers and sources
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
  # Use the pkg_info.py script to determine the module output
  # directory (typically PyTrilinos, but can be PyTrilinos/NOX, etc.)
  EXECUTE_PROCESS(COMMAND
    ${CMAKE_CURRENT_BINARY_DIR}/../util/pkg_info.py -o ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}.i
    OUTPUT_VARIABLE OUTPUT_DIRECTORY
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  SET(CMAKE_SWIG_OUTDIR ${OUTPUT_DIRECTORY})
  # Alternate versions of the module name
  STRING(TOUPPER ${MODULE} MODULE_UPPER)
  STRING(REPLACE "." "_" MODULE_NODOT ${MODULE})
  EXECUTE_PROCESS(COMMAND
    ${CMAKE_CURRENT_BINARY_DIR}/../util/pkg_info.py -m ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}.i
    OUTPUT_VARIABLE MODULE_NAME
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  # The SWIG interface file is for C++
  SET_SOURCE_FILES_PROPERTIES(${MODULE}.i PROPERTIES
    CPLUSPLUS ON
    )
  # Add the SWIG module
  SWIG_ADD_MODULE(${MODULE_NAME} python ${MODULE}.i)
  SWIG_LINK_LIBRARIES(${MODULE_NAME} pytrilinos)
  SWIG_LINK_LIBRARIES(${MODULE_NAME} ${${MODULE_UPPER}_LIBRARIES})
  SET_TARGET_PROPERTIES(${SWIG_MODULE_${MODULE_NAME}_REAL_NAME} PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${OUTPUT_DIRECTORY}
    )
  SET_TARGET_PROPERTIES(${SWIG_MODULE_${MODULE_NAME}_REAL_NAME} PROPERTIES
    LINK_FLAGS "${EXTRA_LINK_ARGS}"
    )
ENDFOREACH(MODULE ${PyTrilinos_MODULES})

#
# Add the PyTrilinos subdirectory, which is top-level package
# directory
ADD_SUBDIRECTORY(PyTrilinos)
