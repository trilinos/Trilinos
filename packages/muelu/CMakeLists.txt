#
# A) Define the package
#

TRIBITS_PACKAGE(MueLu ENABLE_SHADOWING_WARNINGS)

#
# B) Set up package-specific options
#

ASSERT_DEFINED(Xpetra_ENABLE_Epetra)
ASSERT_DEFINED(Xpetra_ENABLE_Tpetra)
ASSERT_DEFINED(Xpetra_INT_LONG_LONG)
ASSERT_DEFINED(Xpetra_ENABLE_Experimental)

ASSERT_DEFINED(Tpetra_INST_DOUBLE)
ASSERT_DEFINED(Tpetra_INST_INT_INT)
ASSERT_DEFINED(Tpetra_INST_INT_LONG)
ASSERT_DEFINED(Tpetra_INST_INT_LONG_LONG)
ASSERT_DEFINED(Tpetra_INST_COMPLEX_DOUBLE)
ASSERT_DEFINED(Tpetra_INST_FLOAT)
ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)

# You can only get Avatar one way
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_AvatarT)
IF( TPL_ENABLE_AVATAR AND ${PACKAGE_NAME}_ENABLE_AvatarT)
  MESSAGE(FATAL_ERROR "MueLu does not support Avatar being built both as a TPL and as an external package.")
ENDIF()

IF ((CMAKE_CXX_COMPILER_ID MATCHES "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "4.7") AND NOT TPL_ENABLE_CUDA)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
ENDIF()

# The following Tribits macro sets MueLu_ENABLE_EXPLICIT_INSTANTIATION
# to be true if ${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION is true.
TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()
TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_DEPRECATED_CODE
  HAVE_${PACKAGE_NAME_UC}_DEPRECATED_CODE
  "Whether MueLu enables deprecated code (that is, anything marked with the MUELU_DEPRECATED macro) at compile time.  Default is ON (deprecated code enabled).  If OFF, then deprecated code does not exist."
  YES
)
IF (HAVE_${PACKAGE_NAME_UC}_DEPRECATED_CODE)
  MESSAGE(STATUS "MueLu: Enabling deprecated code")
ELSE ()
  MESSAGE(STATUS "MueLu: Disabling deprecated code")
ENDIF ()

IF (NOT ${PACKAGE_NAME}_ENABLE_DEPRECATED_CODE AND ${PACKAGE_NAME}_ENABLE_DEPRECATED_TESTS)
  MESSAGE(FATAL_ERROR "MueLu: You cannot enable deprecated tests when deprecated codes is disabled!")
ENDIF ()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_DEPRECATED_TESTS
  HAVE_${PACKAGE_NAME_UC}_DEPRECATED_TESTS
  "Whether MueLu compiles and run tests (or areas of a test) that contain deprecated code (that is, anything marked with the MUELU_DEPRECATED macro).  Default is OFF (deprecated tests disabled).  If ON, then deprecated tests are compiled and tested."
  NO
)

IF (HAVE_${PACKAGE_NAME_UC}_DEPRECATED_TESTS)
  MESSAGE(STATUS "MueLu: Enabling deprecated tests")
ENDIF ()

IF (${PACKAGE_NAME}_ENABLE_Epetra AND NOT ${PACKAGE_NAME}_ENABLE_EpetraExt)
  MESSAGE(FATAL_ERROR "You have enabled Epetra, but not EpetraExt. MueLu requires that either both are enabled, or both are disabled. Please either disable Epetra, or enable EpetraExt.")
ENDIF()


# If you want MueLu Epetra, you need Xpetra Epetra
ASSERT_DEFINED(Xpetra_ENABLE_Epetra)
IF (NOT Xpetra_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Epetra)
  MESSAGE(FATAL_ERROR "MueLu Epetra support requires Xpetra Epetra support")
ENDIF()

# If you want MueLu Epetra, you need Xpetra Tpetra
ASSERT_DEFINED(Xpetra_ENABLE_Tpetra)
IF (NOT Xpetra_ENABLE_Tpetra)
  MESSAGE(FATAL_ERROR "MueLu Tpetra support requires Xpetra Tpetra support")
ENDIF()

# Need to do this because Xpetra_ENABLE_Epetra might have been disabled after
# HAVE_MUELU_EPETRA was first set automatically
GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_EPETRA ${${PACKAGE_NAME}_ENABLE_Epetra})

# Print a warning when we might be lacking a coarse solver.
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos)
IF (${PACKAGE_NAME}_ENABLE_Epetra AND NOT ${PACKAGE_NAME}_ENABLE_Amesos)
  MESSAGE(WARNING "MueLu's Epetra backend is enabled, but Amesos is disabled. This means only iterative coarse solves are available when using Epetra.")
ENDIF ()
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos2)
IF (NOT ${PACKAGE_NAME}_ENABLE_Amesos2)
  MESSAGE(WARNING "MueLu's Amesos2 backend is disabled. This means only iterative coarse solves are available when using Tpetra.")
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
  MESSAGE(STATUS "MueLu: Enabling experimental code")
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

IF (DEFINED ${PACKAGE_NAME}_ENABLE_Kokkos_Refactor)
  MESSAGE(WARNING "MueLu no longer uses the option \'MueLu_ENABLE_Kokkos_Refactor\' and will always build with Kokkos code enabled. \'MueLu_ENABLE_Kokkos_Refactor\' can safely be removed from the CMake invocation.")
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

GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_DEFAULT_GO_LONG OFF)
GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_DEFAULT_GO_LONGLONG OFF)

IF(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
  GLOBAL_SET(${PACKAGE_NAME}_HAVE_GO_INT OFF)
  GLOBAL_SET(${PACKAGE_NAME}_HAVE_GO_LONG OFF)
  GLOBAL_SET(${PACKAGE_NAME}_HAVE_GO_LONG_LONG OFF)

  # <double, int, int>
  # Examples and tests need at least this one, also if Epetra is enabled
  IF(Tpetra_INST_DOUBLE AND Tpetra_INST_INT_INT)
    GLOBAL_SET(${PACKAGE_NAME}_HAVE_GO_INT ON)
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
    GLOBAL_SET(${PACKAGE_NAME}_HAVE_GO_LONG ON)
    GLOBAL_SET(${PACKAGE_NAME}_INST_DOUBLE_INT_LONGINT ON)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_DOUBLE_INT_LONGINT ON)
  ELSE()
    GLOBAL_SET(${PACKAGE_NAME}_INST_DOUBLE_INT_LONGINT OFF)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_DOUBLE_INT_LONGINT OFF)
  ENDIF()


  # <double, int, long long int>
  IF(Tpetra_INST_DOUBLE AND Tpetra_INST_INT_LONG_LONG)
    GLOBAL_SET(${PACKAGE_NAME}_HAVE_GO_LONG_LONG ON)
    GLOBAL_SET(${PACKAGE_NAME}_INST_DOUBLE_INT_LONGLONGINT ON)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_DOUBLE_INT_LONGLONGINT ON)
  ELSE()
    GLOBAL_SET(${PACKAGE_NAME}_INST_DOUBLE_INT_LONGLONGINT OFF)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_DOUBLE_INT_LONGLONGINT OFF)
  ENDIF()

  # <complex, int, int>
  IF(Tpetra_INST_COMPLEX_DOUBLE AND Tpetra_INST_INT_INT)
    GLOBAL_SET(${PACKAGE_NAME}_HAVE_GO_INT ON)
    GLOBAL_SET(${PACKAGE_NAME}_INST_COMPLEX_INT_INT ON)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_COMPLEX_INT_INT ON)
  ELSE()
    GLOBAL_SET(${PACKAGE_NAME}_INST_COMPLEX_INT_INT OFF)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_COMPLEX_INT_INT OFF)
  ENDIF()

  # <complex, int, long long int>
  IF(Tpetra_INST_COMPLEX_DOUBLE AND Tpetra_INST_INT_LONG_LONG)
    GLOBAL_SET(${PACKAGE_NAME}_HAVE_GO_LONG_LONG ON)
    GLOBAL_SET(${PACKAGE_NAME}_INST_COMPLEX_INT_LONG_LONG ON)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_COMPLEX_INT_LONG_LONG ON)
  ELSE()
    GLOBAL_SET(${PACKAGE_NAME}_INST_COMPLEX_INT_LONG_LONG OFF)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_COMPLEX_INT_LONG_LONG OFF)
  ENDIF()

  # <float, int, int>
  IF(Tpetra_INST_FLOAT AND Tpetra_INST_INT_INT)
    GLOBAL_SET(${PACKAGE_NAME}_HAVE_GO_INT ON)
    GLOBAL_SET(${PACKAGE_NAME}_INST_FLOAT_INT_INT ON)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_FLOAT_INT_INT ON)
  ELSE()
    GLOBAL_SET(${PACKAGE_NAME}_INST_FLOAT_INT_INT OFF)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_FLOAT_INT_INT OFF)
  ENDIF()

# <float, int, long long>
  IF(Tpetra_INST_FLOAT AND Tpetra_INST_INT_LONG_LONG)
    GLOBAL_SET(${PACKAGE_NAME}_HAVE_GO_LONG_LONG ON)
    GLOBAL_SET(${PACKAGE_NAME}_INST_FLOAT_INT_LONG_LONG ON)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_FLOAT_INT_LONG_LONG ON)
  ELSE()
    GLOBAL_SET(${PACKAGE_NAME}_INST_FLOAT_INT_LONG_LONG OFF)
    GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_FLOAT_INT_LONG_LONG OFF)
  ENDIF()

  MESSAGE(STATUS "${PACKAGE_NAME}: Enabling ETI support")
  MESSAGE(STATUS "   <float,   int, int>       : ${${PACKAGE_NAME}_INST_FLOAT_INT_INT}")
  MESSAGE(STATUS "   <float,   int, long long> : ${${PACKAGE_NAME}_INST_FLOAT_INT_LONG_LONG}")
  MESSAGE(STATUS "   <double,  int, int>       : ${${PACKAGE_NAME}_INST_DOUBLE_INT_INT}")
  MESSAGE(STATUS "   <double,  int, long>      : ${${PACKAGE_NAME}_INST_DOUBLE_INT_LONGINT}")
  MESSAGE(STATUS "   <double,  int, long long> : ${${PACKAGE_NAME}_INST_DOUBLE_INT_LONGLONGINT}")
  MESSAGE(STATUS "   <complex, int, int>       : ${${PACKAGE_NAME}_INST_COMPLEX_INT_INT}")
  MESSAGE(STATUS "   <complex, int, long long> : ${${PACKAGE_NAME}_INST_COMPLEX_INT_LONG_LONG}")

  # GO preference: int > long long > long
  IF(NOT ${${PACKAGE_NAME}_HAVE_GO_INT})
    # No GO=int (first choice for default)
    IF(NOT ${${PACKAGE_NAME}_HAVE_GO_LONG_LONG})
      #No GO=long long (second choice), must use long
      GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_DEFAULT_GO_LONG ON)
      MESSAGE(STATUS "MueLu: Default GO: long")
    ELSE()
      GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_DEFAULT_GO_LONGLONG ON)
      MESSAGE(STATUS "MueLu: Default GO: long long")
    ENDIF()
  ELSE()
      MESSAGE(STATUS "MueLu: Default GO: int")
  ENDIF()

  # Determine which Tpetra Nodes are used
  GLOBAL_SET (HAVE_${PACKAGE_NAME_UC}_SERIAL  ${Tpetra_INST_SERIAL})
  GLOBAL_SET (HAVE_${PACKAGE_NAME_UC}_OPENMP  ${Tpetra_INST_OPENMP})
  GLOBAL_SET (HAVE_${PACKAGE_NAME_UC}_CUDA    ${Tpetra_INST_CUDA})
  GLOBAL_SET (HAVE_${PACKAGE_NAME_UC}_HIP     ${Tpetra_INST_HIP})
  GLOBAL_SET (HAVE_${PACKAGE_NAME_UC}_SYCL    ${Tpetra_INST_SYCL})
  IF (${PACKAGE_NAME}_ENABLE_Epetra)
    # If Epetra is active, always activate SerialNode
    GLOBAL_SET (HAVE_${PACKAGE_NAME_UC}_SERIAL ON)
  ENDIF ()
  MESSAGE(STATUS "   HAVE_${PACKAGE_NAME_UC}_SERIAL       : ${HAVE_${PACKAGE_NAME_UC}_SERIAL}")
  MESSAGE(STATUS "   HAVE_${PACKAGE_NAME_UC}_OPENMP       : ${HAVE_${PACKAGE_NAME_UC}_OPENMP}")
  MESSAGE(STATUS "   HAVE_${PACKAGE_NAME_UC}_CUDA         : ${HAVE_${PACKAGE_NAME_UC}_CUDA}")
  MESSAGE(STATUS "   HAVE_${PACKAGE_NAME_UC}_HIP          : ${HAVE_${PACKAGE_NAME_UC}_HIP}")
  MESSAGE(STATUS "   HAVE_${PACKAGE_NAME_UC}_SYCL         : ${HAVE_${PACKAGE_NAME_UC}_SYCL}")
ELSE()
  # ETI is turned off
  GLOBAL_SET(${PACKAGE_NAME}_INST_DOUBLE_INT_LONGLONGINT OFF)
  GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_DOUBLE_INT_LONGLONGINT OFF)
  GLOBAL_SET(${PACKAGE_NAME}_INST_COMPLEX_INT_INT OFF)
  GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_COMPLEX_INT_INT OFF)
  GLOBAL_SET(${PACKAGE_NAME}_INST_COMPLEX_INT_LONG_LONG OFF)
  GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_COMPLEX_INT_LONG_LONG OFF)
  GLOBAL_SET(${PACKAGE_NAME}_INST_FLOAT_INT_INT OFF)
  GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_FLOAT_INT_INT OFF)
  GLOBAL_SET(${PACKAGE_NAME}_INST_FLOAT_INT_LONG_LONG OFF)
  GLOBAL_SET(HAVE_${PACKAGE_NAME_UC}_INST_FLOAT_INT_LONG_LONG OFF)

  # If Epetra only is enabled, ETI is off by default. We have to make sure that
  # the HAVE_MUELU_SERIAL variable is set such that the IfpackSmoother and other
  # purely Epetra-related routines are enabled properly
  IF (${PACKAGE_NAME}_ENABLE_Epetra)
    # If Epetra only is active, only activate SerialNode by default
    GLOBAL_SET (HAVE_${PACKAGE_NAME_UC}_SERIAL  ON)
  ENDIF()

  # With ETI off, free to leave int as the default GO
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
ASSERT_DEFINED(Tpetra_INST_COMPLEX_DOUBLE)
IF(NOT Tpetra_INST_COMPLEX_DOUBLE)
  SET(HAVE_${PACKAGE_NAME_UC}_COMPLEX OFF)
ENDIF()

# The HAVE_MUELU_TPETRA_INST_INT_INT flag is enabled if Tpetra is enabled
# with TPETRA_INST_INT_INT. This flag is used to control whether we need
# special stub implementations for the MueLu::Tpetra stack with GO=int.
# taw 9 Oct 2015: This is independent from ETI on/off
ASSERT_DEFINED(Tpetra_INST_INT_INT)
SET (${PACKAGE_NAME}_TPETRA_INST_INT_INT ${Tpetra_INST_INT_INT})
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
  IF(NOT Tpetra_INST_INT_INT OR Tpetra_INST_INT_LONG OR Tpetra_INST_INT_LONG_LONG)
    MESSAGE(FATAL_ERROR "Muemex interfaces require that Tpetra only be instantiated with Tpetra_INST_INT_INT, not Tpetra_INST_INT_LONG or Tpetra_INST_INT_LONG_LONG")
  ENDIF()

  IF (NOT ${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
    MESSAGE(FATAL_ERROR "Muemex interfaces require \"-D${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON\".")
  ENDIF()

  IF(NOT ${PACKAGE_NAME}_ENABLE_Epetra)
    MESSAGE(FATAL_ERROR "Muemex interfaces require \"-D${PACKAGE_NAME}_ENABLE_Epetra:BOOL=ON\".")
  ENDIF()
ENDIF()

if(TPL_ENABLE_MATLAB)
  ADD_SUBDIRECTORY(matlab)
ENDIF()
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
ADD_SUBDIRECTORY(research)
TRIBITS_ADD_TEST_DIRECTORIES(test)

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
#
#
