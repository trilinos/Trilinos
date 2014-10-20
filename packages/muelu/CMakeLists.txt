INCLUDE(CombinedOption)

#
# A) Define the package
#

TRIBITS_PACKAGE(MueLu ENABLE_SHADOWING_WARNINGS)

#
# B) Set up package-specific options
#
IF ((CMAKE_CXX_COMPILER_ID MATCHES "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "4.7"))
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
ENDIF()


TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()
TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

IF (${PACKAGE_NAME}_ENABLE_Epetra AND NOT ${PACKAGE_NAME}_ENABLE_EpetraExt)
  MESSAGE(FATAL_ERROR "You have enabled Epetra, but not EpetraExt. MueLu requires that either both are enabled, or both are disabled. Please either disable Epetra, or enable EpetraExt.")
ENDIF()

# Disable Ifpack2 if Tpetra is off
IF (NOT ${PACKAGE_NAME}_ENABLE_Tpetra)
  SET(${PACKAGE_NAME}_ENABLE_Ifpack2 OFF)
  GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_IFPACK2 OFF)
ENDIF ()

# Enable timer synchronization via barriers
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_TIMER_SYNCHRONIZATION
  HAVE_${PACKAGE_NAME_UC}_TIMER_SYNCHRONIZATION
  "Enable timer synchronization in MueLu via barriers."
  NO)

IF (${PACKAGE_NAME}_ENABLE_TIMER_SYNCHRONIZATION)
  MESSAGE(STATUS "Enabling timer synchronization")
ENDIF ()

# Enable experimental code
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Experimental
  HAVE_${PACKAGE_NAME_UC}_EXPERIMENTAL
  "Enable experimental code."
  NO)
IF (${PACKAGE_NAME}_ENABLE_Experimental AND NOT Xpetra_ENABLE_Experimental)
  MESSAGE(FATAL_ERROR "Option MueLu_ENABLE_Experimental=ON requires Xpetra_ENABLE_Experimental=ON.")
ENDIF ()

IF (${PACKAGE_NAME}_ENABLE_Experimental)
  MESSAGE(STATUS "Enabling MueLu experimental code")
ENDIF()

# Memory testing
IF (${PACKAGE_NAME}_ENABLE_MEMORY_TESTING)
  SET(MUELU_VALGRIND_TESTING YES)
ELSE()
  SET(MUELU_VALGRIND_TESTING NO)
ENDIF()

# Profiling
TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_PROFILING
  HAVE_${PACKAGE_NAME_UC}_PROFILING
  "Enable detailed timing and memory profiling of MueLu."
  OFF)

IF (${PACKAGE_NAME}_ENABLE_PROFILING)
  MESSAGE(STATUS "Enabling MueLu profiling code")
ENDIF()

# Profiling (expert users only. There is no TPL for Google Perftools)
TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_GOOGLE_PERFTOOLS
  HAVE_${PACKAGE_NAME_UC}_GOOGLE_PERFTOOLS
  "Enable detailed memory profiling using Google Perftools."
  OFF)

# Boost
# We don't actually need this option, but boost is not *really* needed, and
# I really want to hide unnecessary warnings
TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Boost_for_real
  HAVE_${PACKAGE_NAME_UC}_BOOST_FOR_REAL
  "Enable Boost"
  OFF)

# Explicit instantiation options:
# - MueLu_ENABLE_EXPLICIT_INSTANTIATION
# - MueLu_INST_DOUBLE_INT_INT         (= ON by default)
# - MueLu_INST_DOUBLE_INT_LONGLONGINT (= ON by default if TEUCHOS_ENABLE_LONG_LONG_INT enabled)
# - MueLu_INST_COMPLEX_INT_INT        (= ON by default if TEUCHOS_ENABLE_COMPLEX enabled)       TODO: we also need to check Tpetra_COMPLEX if EI of Tpetra=ON
IF(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)

  # <double, int, int>
  # Examples and tests need at least this one
  TRIBITS_ADD_OPTION_AND_DEFINE(
    ${PACKAGE_NAME}_INST_DOUBLE_INT_INT
    HAVE_${PACKAGE_NAME_UC}_INST_DOUBLE_INT_INT
    "Enable the use of explicit template instantiation with <double, int, int>."
    ON
    )
  IF(${PACKAGE_NAME}_ENABLE_EXAMPLES OR ${PACKAGE_NAME}_ENABLE_TESTS)
    IF (NOT ${PACKAGE_NAME}_INST_DOUBLE_INT_INT)
      MESSAGE(FATAL_ERROR "Error: MueLu examples and tests needs explicit template instantation with <double, int, int>. Turn on the option MueLu_INST_DOUBLE_INT_INT or disable examples and tests.")
    ENDIF()
  ENDIF()

  # <double, int, long long int>
  IF(Teuchos_ENABLE_LONG_LONG_INT)
    SET(${PACKAGE_NAME}_INST_DOUBLE_INT_LONGLONGINT_DEFAULT ON)
  ELSE()
    SET(${PACKAGE_NAME}_INST_DOUBLE_INT_LONGLONGINT_DEFAULT OFF)
  ENDIF()
  TRIBITS_ADD_OPTION_AND_DEFINE(
    ${PACKAGE_NAME}_INST_DOUBLE_INT_LONGLONGINT
    HAVE_${PACKAGE_NAME_UC}_INST_DOUBLE_INT_LONGLONGINT
    "Enable the use of explicit template instantiation with <double, int, long long int>."
    ${${PACKAGE_NAME}_INST_DOUBLE_INT_LONGLONGINT_DEFAULT}
    )

  # <complex, int, int>
  # Teuchos_ENABLE_COMPLEX is ON by default (in release mode or when SS is enabled) but we don't want MueLu to be instantiated explicitly with complex by default as it is not frequently used.
  #  IF(Teuchos_ENABLE_COMPLEX)
  #    SET(${PACKAGE_NAME}_INST_COMPLEX_INT_INT_DEFAULT ON)
  #  ELSE()
    SET(${PACKAGE_NAME}_INST_COMPLEX_INT_INT_DEFAULT OFF)
  #  ENDIF()
  TRIBITS_ADD_OPTION_AND_DEFINE(
    ${PACKAGE_NAME}_INST_COMPLEX_INT_INT
    HAVE_${PACKAGE_NAME_UC}_INST_COMPLEX_INT_INT
    "Enable the use of explicit template instantiation with <complex, int, int>."
    ${${PACKAGE_NAME}_INST_COMPLEX_INT_INT_DEFAULT}
    )
  MESSAGE(STATUS "${PACKAGE_NAME} ETI support enabled")
  MESSAGE(STATUS "<double,  int, int>       : ${${PACKAGE_NAME}_INST_DOUBLE_INT_INT}")
  MESSAGE(STATUS "<double,  int, long long> : ${${PACKAGE_NAME}_INST_DOUBLE_INT_LONGLONGINT}")
  MESSAGE(STATUS "<complex, int, int>       : ${${PACKAGE_NAME}_INST_COMPLEX_INT_INT}")

ENDIF()

# HAVE_MUELU_COMPLEX=ON iff:
# - Teuchos_ENABLE_COMPLEX=ON and Explicit Instantiation = OFF
# or
# - HAVE_MUELU_INST_COMPLEX_INT_INT
#
# TODO: do we want also a MueLu_ENABLE_Complex option that turns ON automatically what we need?
IF(((NOT ${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION) AND Teuchos_ENABLE_COMPLEX)
    OR
    (HAVE_MUELU_INST_COMPLEX_INT_INT))
  SET(HAVE_${PACKAGE_NAME_UC}_COMPLEX ON)
ELSE()
  SET(HAVE_${PACKAGE_NAME_UC}_COMPLEX OFF)
ENDIF()


# ML's MMM
TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_ML_MMM
  HAVE_${PACKAGE_NAME_UC}_ML_MMM
  "Use ML's Matrix-Matrix-Multiply when available"
  OFF)
IF (${PACKAGE_NAME}_ENABLE_ML_MMM AND NOT ${PACKAGE_NAME}_ENABLE_ML)
  MESSAGE(FATAL_ERROR "Option MueLu_ENABLE_ML_MMM=ON requires MueLu_ENABLE_ML=ON.")
ENDIF ()


#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(adapters)
#address bugzilla bug#6236
IF (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
    AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/doc/CMakeLists.txt)
  ADD_SUBDIRECTORY(doc)
ENDIF()

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(research)
TRIBITS_ADD_TEST_DIRECTORIES(test)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
#
#
