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
# - Enabling/disabling of scalar/ordinal types is controlled by Tpetra
IF(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)

  # <double, int, int>
  # Examples and tests need at least this one, also if Epetra is enabled
  IF(Tpetra_INST_DOUBLE AND Tpetra_INST_INT_INT)
    GLOBAL_SET(${PACKAGE_NAME}_INST_DOUBLE_INT_INT ON)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_DOUBLE_INT_INT ON)
  ELSE()
    GLOBAL_SET(${PACKAGE_NAME}_INST_DOUBLE_INT_INT OFF)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_DOUBLE_INT_INT OFF)
  ENDIF()

  IF(${PACKAGE_NAME}_ENABLE_EXAMPLES OR ${PACKAGE_NAME}_ENABLE_TESTS)
    IF (NOT ${PACKAGE_NAME}_INST_DOUBLE_INT_INT)
      MESSAGE(FATAL_ERROR "Error: MueLu examples and tests needs explicit template instantation with <double, int, int>. Turn on the options Tpetra_INST_DOUBLE and Tpetra_INST_INT_INT or disable examples and tests.")
    ENDIF()
  ENDIF()

  IF(${PACKAGE_NAME}_ENABLE_Epetra)
    IF (NOT ${PACKAGE_NAME}_INST_DOUBLE_INT_INT)
      MESSAGE(FATAL_ERROR "Error: MueLu needs explicit template instantation with <double, int, int> when enabling Epetra. Turn on the options Tpetra_INST_DOUBLE and Tpetra_INST_INT_INT or disable Epetra.")
    ENDIF()
  ENDIF()

  # <double, int, long>
  IF(Tpetra_INST_DOUBLE AND Tpetra_INST_INT_LONG)
    GLOBAL_SET(${PACKAGE_NAME}_INST_DOUBLE_INT_LONGINT ON)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_DOUBLE_INT_LONGINT ON)
  ELSE()
    GLOBAL_SET(${PACKAGE_NAME}_INST_DOUBLE_INT_LONGINT OFF)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_DOUBLE_INT_LONGINT OFF)
  ENDIF()


  # <double, int, long long int>
   IF(Tpetra_INST_DOUBLE AND Tpetra_INST_INT_LONG_LONG)
    GLOBAL_SET(${PACKAGE_NAME}_INST_DOUBLE_INT_LONGLONGINT ON)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_DOUBLE_INT_LONGLONGINT ON)
  ELSE()
    GLOBAL_SET(${PACKAGE_NAME}_INST_DOUBLE_INT_LONGLONGINT OFF)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_DOUBLE_INT_LONGLONGINT OFF)
  ENDIF()

  # <complex, int, int>
  IF(Tpetra_INST_COMPLEX_DOUBLE AND Tpetra_INST_INT_INT)
    GLOBAL_SET(${PACKAGE_NAME}_INST_COMPLEX_INT_INT ON)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_COMPLEX_INT_INT ON)
  ELSE()
    GLOBAL_SET(${PACKAGE_NAME}_INST_COMPLEX_INT_INT OFF)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_COMPLEX_INT_INT OFF)
  ENDIF()

  # <float, int, int>
  IF(Tpetra_INST_FLOAT AND Tpetra_INST_INT_INT)
    GLOBAL_SET(${PACKAGE_NAME}_INST_FLOAT_INT_INT ON)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_FLOAT_INT_INT ON)
  ELSE()
    GLOBAL_SET(${PACKAGE_NAME}_INST_FLOAT_INT_INT OFF)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_FLOAT_INT_INT OFF)
  ENDIF()

  MESSAGE(STATUS "${PACKAGE_NAME} ETI support enabled")
  MESSAGE(STATUS "<float,   int, int>       : ${${PACKAGE_NAME}_INST_FLOAT_INT_INT}")
  MESSAGE(STATUS "<double,  int, int>       : ${${PACKAGE_NAME}_INST_DOUBLE_INT_INT}")
  MESSAGE(STATUS "<double,  int, long>      : ${${PACKAGE_NAME}_INST_DOUBLE_INT_LONGINT}")
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
