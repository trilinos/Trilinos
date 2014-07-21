# -*- cmake -*-
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# Define the package
TRIBITS_PACKAGE(Domi DISABLE_STRONG_WARNINGS)

#
# Standardized debugging options
TRIBITS_ADD_DEBUG_OPTION()

#
# Set up package-specific options
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_ABC
  HAVE_DOMI_ARRAY_BOUNDSCHECK
  "Enable Domi feature for array bounds checking."
  ${${PROJECT_NAME}_ENABLE_DEBUG}
  )

########################################################################

# Python interpreter, library, and header files are required.
#SET(CMAKE_MODULE_PATH
#  ${CMAKE_PROJECT_DIR}/packages/PyTrilinos/cmake/FindPythonInclude.cmake
#  ${CMAKE_MODULE_PATH})
SET(CMAKE_PREFIX_PATH
  ${CMAKE_PROJECT_DIR}/packages/PyTrilinos/cmake
  ${CMAKE_PREFIX_PATH})
FIND_PACKAGE(PythonInterp 2.2 REQUIRED)
FIND_PACKAGE(PythonLibs 2.2 REQUIRED)
FIND_PACKAGE(PythonInclude REQUIRED)

# Python module numpy is required
FIND_PACKAGE(NumPy 1.0 REQUIRED)

# Simple wrapper interface generator (SWIG) is required
SET(SWIG_VERSION_FLOOR   2.0.0)
SET(SWIG_VERSION_CEILING 2.0.8)
FIND_PACKAGE(SWIG REQUIRED)
IF(${SWIG_VERSION} VERSION_LESS ${SWIG_VERSION_FLOOR})
  MESSAGE(FATAL_ERROR
    "SWIG version " ${SWIG_VERSION}
    " is less than minimum version " ${SWIG_VERSION_FLOOR}
    )
ENDIF()
IF(${SWIG_VERSION} VERSION_GREATER ${SWIG_VERSION_CEILING})
  MESSAGE(FATAL_ERROR
    "SWIG version " ${SWIG_VERSION}
    " is greater than maximum version " ${SWIG_VERSION_CEILING}
    )
ENDIF()
# ${SWIG_USE_FILE} is hard-coded to be the installed UseSWIG.cmake
# file.  We want to include the local, modified version.
INCLUDE(UseSWIG)

########################################################################

#
# Add the libraries, tests, and examples
ADD_SUBDIRECTORY(src)
TRIBITS_ADD_TEST_DIRECTORIES(test)

#
# Add the Parameter List documentation
TRIBITS_ADD_TEST_DIRECTORIES(doc/parameterList)

#
# Do standard postprocessing
TRIBITS_PACKAGE_POSTPROCESS()
