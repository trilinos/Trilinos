INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)
INCLUDE(CombinedOption)

#
# A) Define the package
#

PACKAGE(MueLu)

#
# B) Set up package-specific options
#

PACKAGE_ADD_EXPLICIT_INSTANTIATION_OPTION()
PACKAGE_ADD_DEBUG_OPTION()
PACKAGE_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

# Profiling
ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_PROFILING
  HAVE_${PACKAGE_NAME_UC}_PROFILING
  "Enable detailed timing and memory profiling of MueLu."
  OFF)

# Explicit instantiation options:
# - MueLu_ENABLE_EXPLICIT_INSTANTIATION
# - MueLu_ENABLE_EI_DOUBLE_INT_INT         (= ON by default)
# - MueLu_ENABLE_EI_DOUBLE_INT_LONGLONGINT (= ON by default if TEUCHOS_ENABLE_LONG_LONG_INT enabled)
IF(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)

  # <double, int, int>
  # Examples and tests needs at least this one
  ADD_OPTION_AND_DEFINE(
    ${PACKAGE_NAME}_ENABLE_EI_DOUBLE_INT_INT
    HAVE_${PACKAGE_NAME_UC}_EI_DOUBLE_INT_INT
    "Enable the use of explicit template instantiation with <double, int, int>."
    ON
    )
  IF(${PACKAGE_NAME}_ENABLE_EXAMPLES OR ${PACKAGE_NAME}_ENABLE_TESTS) 
    IF (NOT ${PACKAGE_NAME}_ENABLE_EI_DOUBLE_INT_INT)
      MESSAGE(FATAL_ERROR "Error: MueLu examples and tests needs explicit template instantation with <double, int, int>. Turn on the option MueLu_ENABLE_EI_DOUBLE_INT_INT or disable examples and tests.")
    ENDIF()
  ENDIF()

  # <double, int, long long int>
  IF(${Teuchos_ENABLE_LONG_LONG_INT})
    SET(${PACKAGE_NAME}_ENABLE_EI_DOUBLE_INT_LONGLONGINT_DEFAULT ON)
  ELSE()
    SET(${PACKAGE_NAME}_ENABLE_EI_DOUBLE_INT_LONGLONGINT_DEFAULT OFF)
  ENDIF()
  ADD_OPTION_AND_DEFINE(
    ${PACKAGE_NAME}_ENABLE_EI_DOUBLE_INT_LONGLONGINT
    HAVE_${PACKAGE_NAME_UC}_EI_DOUBLE_INT_LONGLONGINT
    "Enable the use of explicit template instantiation with <double, int, long long int>."
    ${PACKAGE_NAME}_ENABLE_EI_DOUBLE_INT_LONGLONGINT_DEFAULT
    )
ENDIF()

#
# Add directory xpetra/cmake/ to the cmake module search path
#

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PACKAGE_SOURCE_DIR}/xpetra/cmake)

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(gallery) # iff EXAMPLE or TEST ENABLED
ADD_SUBDIRECTORY(xpetra)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(adapters)
ADD_SUBDIRECTORY(doc)
ADD_SUBDIRECTORY(utils)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)
PACKAGE_ADD_TEST_DIRECTORIES(test)

#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
