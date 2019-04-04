INCLUDE(CombinedOption)

#
# A) Define the package
#

TRIBITS_PACKAGE(MueLu ENABLE_SHADOWING_WARNINGS)

#
# B) Set up package-specific options
#

IF (NOT ${PACKAGE_NAME}_ENABLE_Epetra AND NOT ${PACKAGE_NAME}_ENABLE_Tpetra)
  MESSAGE(FATAL_ERROR "Neither Epetra nor Tpetra are enabled.  MueLu requires at least one of them to be enabled.")
ENDIF()

ASSERT_DEFINED(Xpetra_INT_LONG_LONG)
ASSERT_DEFINED(Xpetra_ENABLE_Experimental)
ASSERT_DEFINED(Xpetra_ENABLE_Kokkos_Refactor)
IF (${PACKAGE_NAME}_ENABLE_KokkosCore AND NOT ${PACKAGE_NAME}_ENABLE_Tpetra)
  MESSAGE(FATAL_ERROR "MueLu requires that if Kokkos is enabled, Tpetra must also be enabled.")
ENDIF()

if (${PACKAGE_NAME}_ENABLE_Tpetra)
  ASSERT_DEFINED(Tpetra_INST_DOUBLE)
  ASSERT_DEFINED(Tpetra_INST_INT_INT)
  ASSERT_DEFINED(Tpetra_INST_INT_LONG)
  ASSERT_DEFINED(Tpetra_INST_INT_LONG_LONG)
  ASSERT_DEFINED(Tpetra_INST_COMPLEX_DOUBLE)
  ASSERT_DEFINED(Tpetra_INST_FLOAT)
  ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
ENDIF()

# JJH: pulled this logic from Tpetra.  Comments below about CUDA are from Tpetra developers.
# This is redundant, since MueLu also requires Kokkos, and Kokkos requires C++11.
# "I left a special case for CUDA because I'm not sure whether TriBITS
# counts C++11 as enabled in that case; please feel free to revise if
# you know better."
ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_CXX11)
ASSERT_DEFINED(TPL_ENABLE_CUDA)
IF (NOT ${PROJECT_NAME}_ENABLE_CXX11 AND NOT TPL_ENABLE_CUDA)
  MESSAGE(FATAL_ERROR "MueLu requires that C++11 support in Trilinos be enabled.  This should happen automatically if your compiler is C++11 compliant.")
ENDIF ()

IF ((CMAKE_CXX_COMPILER_ID MATCHES "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "4.7") AND NOT TPL_ENABLE_CUDA)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
ENDIF()

# The following Tribits macro sets MueLu_ENABLE_EXPLICIT_INSTANTIATION
# to be true if ${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION is true.
TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()
IF (${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION AND NOT ${PACKAGE_NAME}_ENABLE_Tpetra)
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION  OFF)
  GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_EXPLICIT_INSTANTIATION OFF)
  MESSAGE(STATUS "Disabling MueLu explicit template instantation (ETI) because Tpetra is disabled.")
ENDIF()
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

# Enable additive variant
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_ADDITIVE_VARIANT
  HAVE_${PACKAGE_NAME_UC}_ADDITIVE_VARIANT
  "Enable Additive Variant."
  NO)

# Enable matrix splitting
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_REGION_SPLITTING
  HAVE_${PACKAGE_NAME_UC}_REGION_SPLITTING
  "Enable Region Splitting."
  NO)

# mfh 12 Sep 2017: See #1726.  If Tpetra is enabled, then MueLu's
# Tutorial should not build unless GO=int and Node=OpenMP or Serial
# are all enabled.

SET(${PACKAGE_NAME}_ENABLE_Tutorial_DEFAULT YES)
IF (DEFINED ${PACKAGE_NAME}_ENABLE_Tpetra)
  IF (${PACKAGE_NAME}_ENABLE_Tpetra)
    IF (NOT DEFINED Tpetra_INST_INT_INT OR NOT Tpetra_INST_INT_INT)
      SET(${PACKAGE_NAME}_ENABLE_Tutorial_DEFAULT NO)
    ENDIF ()
    IF (NOT DEFINED Tpetra_INST_SERIAL OR NOT Tpetra_INST_SERIAL)
      IF (NOT DEFINED Tpetra_INST_OPENMP OR NOT Tpetra_INST_OPENMP)
        SET(${PACKAGE_NAME}_ENABLE_Tutorial_DEFAULT NO)
      ENDIF ()
    ENDIF ()
    IF (NOT DEFINED Tpetra_INST_OPENMP OR NOT Tpetra_INST_OPENMP)
      IF (NOT DEFINED Tpetra_INST_SERIAL OR NOT Tpetra_INST_SERIAL)
        MESSAGE (STATUS "MueLu's tutorial requires that either the Serial or OpenMP Nodes in Tpetra be enabled.  Thus, I am disabling MueLu's tutorial by default.")
        SET(${PACKAGE_NAME}_ENABLE_Tutorial_DEFAULT NO)
      ENDIF ()
    ENDIF ()
  ENDIF ()
ENDIF ()

IF (NOT ${PACKAGE_NAME}_ENABLE_Tutorial_DEFAULT)
  MESSAGE (STATUS "MueLu's tutorial requires for a correct build that if Tpetra is enabled, then GlobalOrdinal=int must be enabled, and either the Serial or OpenMP Nodes in Tpetra be enabled.  Since this is not the case, I am disabling MueLu's tutorial by default.  See Trilinos GitHub Issue #1726 for details: https://github.com/trilinos/Trilinos/issues/1726")
ENDIF ()

# cmake option to enable building of MueLu Tutorial
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Tutorial
  ""
  "Enable MueLu tutorial."
  ${PACKAGE_NAME}_ENABLE_Tutorial_DEFAULT)

IF (NOT ${PACKAGE_NAME}_ENABLE_Tutorial_DEFAULT AND ${PACKAGE_NAME}_ENABLE_Tutorial)
  MESSAGE (STATUS "You chose to enable MueLu's tutorial, even though I told you that it may not build correctly.  You are responsible for whatever build or link errors may occur as a result.")
ENDIF ()

# Split ETI cpp files
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_SPLIT_ETI_CPP_FILES
  HAVE_${PACKAGE_NAME_UC}_SPLIT_ETI_CPP_FILES
  "Split ETI cpp files into extra files for each enabled configuration (not recommended)."
  NO)
IF (${PACKAGE_NAME}_ENABLE_SPLIT_ETI_CPP_FILES)
  MESSAGE(STATUS "Split MueLu ETI files: this may slow down the compilation process.")
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

# Define option to enable code that is refactored to use Kokkos
# This is ON by default when Kokkos is enabled.
IF (${PACKAGE_NAME}_ENABLE_KokkosCore)
  SET(Kokkos_Refactor_On_By_Default ON)
  MESSAGE(STATUS "MueLu: Kokkos_Refactor ON by default")
ELSE()
  SET(Kokkos_Refactor_On_By_Default OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Kokkos_Refactor
  HAVE_${PACKAGE_NAME_UC}_KOKKOS_REFACTOR
  "Enable Xpetra code refactored to use Kokkos."
  ${Kokkos_Refactor_On_By_Default})

# Define option to use code that is refactored to use Kokkos by
# default. This is useful for running all non-unit tests with Kokkos
# without having to set "use kokkos refactor" in all parameter lists.
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Kokkos_Refactor_Use_By_Default
  HAVE_${PACKAGE_NAME_UC}_KOKKOS_REFACTOR_USE_BY_DEFAULT
  "Use MueLu code refactored to use Kokkos by default."
  NO)


IF (${PACKAGE_NAME}_ENABLE_Kokkos_Refactor AND NOT Xpetra_ENABLE_Kokkos_Refactor)
  MESSAGE(FATAL_ERROR "Kokkos Refactor in MueLu requires one to also enable Kokkos Refactor in Xpetra (-DXpetra_ENABLE_Kokkos_Refactor:BOOL=ON).")
ENDIF()

# Define HAVE_MUELU_STRATIMIKOS if Stratimikos package is enabled
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Stratimikos)
GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_STRATIMIKOS ${${PACKAGE_NAME}_ENABLE_Stratimikos})

# Define HAVE_MUELU_TEKO if Teko package is enabled
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Teko)
GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_TEKO ${${PACKAGE_NAME}_ENABLE_Teko})

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

  #IF(${PACKAGE_NAME}_ENABLE_EXAMPLES OR ${PACKAGE_NAME}_ENABLE_TESTS)
  #  IF (NOT ${PACKAGE_NAME}_INST_DOUBLE_INT_INT)
  #    MESSAGE(FATAL_ERROR "Error: MueLu examples and tests needs explicit template instantation with <double, int, int>. Turn on the options Tpetra_INST_DOUBLE and Tpetra_INST_INT_INT or disable examples and tests.")
  #  ENDIF()
  #ENDIF()

  #IF(${PACKAGE_NAME}_ENABLE_Epetra)
  #  IF (NOT ${PACKAGE_NAME}_INST_DOUBLE_INT_INT)
  #    MESSAGE(FATAL_ERROR "Error: MueLu needs explicit template instantation with <double, int, int> when enabling Epetra. Turn on the options Tpetra_INST_DOUBLE and Tpetra_INST_INT_INT or disable Epetra.")
  #  ENDIF()
  #ENDIF()

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

  # <complex, int, long long int>
  IF(Tpetra_INST_COMPLEX_DOUBLE AND Tpetra_INST_INT_LONG_LONG)
    GLOBAL_SET(${PACKAGE_NAME}_INST_COMPLEX_INT_LONG_LONG ON)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_COMPLEX_INT_LONG_LONG ON)
  ELSE()
    GLOBAL_SET(${PACKAGE_NAME}_INST_COMPLEX_INT_LONG_LONG OFF)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_COMPLEX_INT_LONG_LONG OFF)
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

  # If Tpetra is enabled, determine which Nodes are used
  IF (${PACKAGE_NAME}_ENABLE_Tpetra)
    GLOBAL_SET (HAVE_${PACKAGE_NAME_UC}_SERIAL  ${Tpetra_INST_SERIAL})
    GLOBAL_SET (HAVE_${PACKAGE_NAME_UC}_PTHREAD ${Tpetra_INST_PTHREAD})
    GLOBAL_SET (HAVE_${PACKAGE_NAME_UC}_OPENMP  ${Tpetra_INST_OPENMP})
    GLOBAL_SET (HAVE_${PACKAGE_NAME_UC}_CUDA    ${Tpetra_INST_CUDA})
  ENDIF()
  IF (${PACKAGE_NAME}_ENABLE_Epetra)
    # If Epetra is active, always activate SerialNode
    GLOBAL_SET (HAVE_${PACKAGE_NAME_UC}_SERIAL ON)
  ENDIF ()
  MESSAGE(STATUS "HAVE_${PACKAGE_NAME_UC}_SERIAL       : ${HAVE_${PACKAGE_NAME_UC}_SERIAL}")
  MESSAGE(STATUS "HAVE_${PACKAGE_NAME_UC}_PTHREAD      : ${HAVE_${PACKAGE_NAME_UC}_PTHREAD}")
  MESSAGE(STATUS "HAVE_${PACKAGE_NAME_UC}_OPENMP       : ${HAVE_${PACKAGE_NAME_UC}_OPENMP}")
  MESSAGE(STATUS "HAVE_${PACKAGE_NAME_UC}_CUDA         : ${HAVE_${PACKAGE_NAME_UC}_CUDA}")
ELSE()
  # ETI is turned off

  # If Epetra only is enabled, ETI is off by default. We have to make sure that
  # the HAVE_MUELU_SERIAL variable is set such that the IfpackSmoother and other
  # purely Epetra-related routines are enabled properly
  IF (${PACKAGE_NAME}_ENABLE_Epetra)
    # If Epetra only is active, only activate SerialNode by default
    GLOBAL_SET (HAVE_${PACKAGE_NAME_UC}_SERIAL  ON)
  ENDIF()
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

# mfh 18 Aug 2015: Revise whether complex is enabled, based on
# Tpetra's settings.  These are independent of whether ETI is enabled!
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Tpetra)
IF(${PACKAGE_NAME}_ENABLE_Tpetra)
  ASSERT_DEFINED(Tpetra_INST_COMPLEX_DOUBLE)
  IF(NOT Tpetra_INST_COMPLEX_DOUBLE)
    SET(HAVE_${PACKAGE_NAME_UC}_COMPLEX OFF)
  ENDIF()
ENDIF()

# The HAVE_MUELU_TPETRA_INST_INT_INT flag is enabled if Tpetra is enabled
# with TPETRA_INST_INT_INT. This flag is used to control whether we need
# special stub implementations for the MueLu::Tpetra stack with GO=int.
# taw 9 Oct 2015: This is independent from ETI on/off
IF (${PACKAGE_NAME}_ENABLE_Tpetra)
  ASSERT_DEFINED(Tpetra_INST_INT_INT)
  SET (${PACKAGE_NAME}_TPETRA_INST_INT_INT ${Tpetra_INST_INT_INT})
ENDIF()
GLOBAL_SET (HAVE_MUELU_TPETRA_INST_INT_INT ${${PACKAGE_NAME}_TPETRA_INST_INT_INT})

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

# The build has to be in a very specific order because of the existence of Matlab interface:
#   1. Build the majority of the code in src/   [except src/Interface; src/CMakeLists.txt is responsible for ignoring that]
#   2. Build matlab factories                   [located in matlab/; matlab/bin is ignored in matlab/CMakeLists.txt]
#   3. Build the interpreter                    [requires all factories]
#   4. Build the adapters                       [required for Matlab]
#   5. Build the rest of matlab code
# This order required splitting of the src/ library into two: muelu and muelu-interface, as steps 1. and 3. are separated.
# A simple diagram explains this:
#                                               /  muelu
# matlab/bin - muelu-adapters - muelu-interface      |
#                                               \  matlab
ADD_SUBDIRECTORY(src)

if(TPL_ENABLE_MATLAB)
  if (NOT ${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
    MESSAGE(FATAL_ERROR "Muemex interfaces require \"-D${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON\".")
  ENDIF()
  IF(NOT ${PACKAGE_NAME}_ENABLE_Epetra OR NOT ${PACKAGE_NAME}_ENABLE_Tpetra)
      MESSAGE(FATAL_ERROR "Muemex interfaces require \"-D${PACKAGE_NAME}_ENABLE_Epetra:BOOL=ON\" and \"-D${PACKAGE_NAME}_ENABLE_Tpetra:BOOL=ON\".")
  ENDIF()
ENDIF()

if(TPL_ENABLE_MATLAB)
  ADD_SUBDIRECTORY(matlab)
ENDIF()
ADD_SUBDIRECTORY(src/Interface)
ADD_SUBDIRECTORY(adapters)
if(TPL_ENABLE_MATLAB)
  ADD_SUBDIRECTORY(matlab/bin)
  ADD_SUBDIRECTORY(matlab/tests)
ENDIF()

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
