# Include the package macros
INCLUDE(PackageMacros)

# Set the package name
PACKAGE(PyTrilinos  DISABLE_STRONG_WARNINGS)

# Set the package version number
SET(${PACKAGE_NAME}_VERSION 4.2d)

# Define the list of PyTrilinos modules
SET(${PACKAGE_NAME}_MODULES Epetra)

IF(PyTrilinos_ENABLE_Teuchos)
  APPEND_SET(${PACKAGE_NAME}_MODULES Teuchos)
ENDIF(PyTrilinos_ENABLE_Teuchos)

IF(PyTrilinos_ENABLE_Triutils)
  APPEND_SET(${PACKAGE_NAME}_MODULES TriUtils)
ENDIF(PyTrilinos_ENABLE_Triutils)

IF(PyTrilinos_ENABLE_EpetraExt)
  APPEND_SET(${PACKAGE_NAME}_MODULES EpetraExt)
ENDIF(PyTrilinos_ENABLE_EpetraExt)

IF(PyTrilinos_ENABLE_AztecOO)
  APPEND_SET(${PACKAGE_NAME}_MODULES AztecOO)
ENDIF(PyTrilinos_ENABLE_AztecOO)

IF(PyTrilinos_ENABLE_Galeri)
  APPEND_SET(${PACKAGE_NAME}_MODULES Galeri)
ENDIF(PyTrilinos_ENABLE_Galeri)

IF(PyTrilinos_ENABLE_Amesos)
  APPEND_SET(${PACKAGE_NAME}_MODULES Amesos)
ENDIF(PyTrilinos_ENABLE_Amesos)

IF(PyTrilinos_ENABLE_Ifpack)
  APPEND_SET(${PACKAGE_NAME}_MODULES IFPACK)
ENDIF(PyTrilinos_ENABLE_Ifpack)

IF(PyTrilinos_ENABLE_Komplex)
  APPEND_SET(${PACKAGE_NAME}_MODULES Komplex)
ENDIF(PyTrilinos_ENABLE_Komplex)

IF(PyTrilinos_ENABLE_Anasazi)
  APPEND_SET(${PACKAGE_NAME}_MODULES Anasazi)
ENDIF(PyTrilinos_ENABLE_Anasazi)

IF(PyTrilinos_ENABLE_Pliris)
  APPEND_SET(${PACKAGE_NAME}_MODULES Pliris)
ENDIF(PyTrilinos_ENABLE_Pliris)

IF(PyTrilinos_ENABLE_ML)
  APPEND_SET(${PACKAGE_NAME}_MODULES ML)
ENDIF(PyTrilinos_ENABLE_ML)

IF(PyTrilinos_ENABLE_NOX)
  APPEND_SET(${PACKAGE_NAME}_MODULES NOX)
ENDIF(PyTrilinos_ENABLE_NOX)

IF(PyTrilinos_ENABLE_LOCA)
  APPEND_SET(${PACKAGE_NAME}_MODULES LOCA)
ENDIF(PyTrilinos_ENABLE_LOCA)

# Print the enabled PyTrilinos packages
MESSAGE(STATUS "Enabled ${PACKAGE_NAME} modules:")
FOREACH(MODULE ${${PACKAGE_NAME}_MODULES})
  MESSAGE(STATUS "    ${MODULE}")
ENDFOREACH(MODULE ${${PACKAGE_NAME}_MODULES})

# Python interpreter and header files are required
FIND_PACKAGE(PythonInterp 2.2 REQUIRED)
FIND_PACKAGE(PythonInclude REQUIRED)

# Python module numpy is required
FIND_PACKAGE(NumPy 1.0 REQUIRED)

# Simple wrapper interface generator (SWIG) is required
FIND_PACKAGE(SWIG 1.3.29 REQUIRED)
INCLUDE(${SWIG_USE_FILE})

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

#PACKAGE_ADD_TEST_DIRECTORIES(test)

#PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

PACKAGE_POSTPROCESS()
