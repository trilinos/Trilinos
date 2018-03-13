#
# Teuchos is both a package and can be configured and built as it own TriBITS
# project!
#

INCLUDE(CombinedOption)


IF (TRIBITS_PROCESSING_PACKAGE)


#
# Processing TriBITS package!
#

#
# A) Forward delcare the package so that certain options are also defined for
# subpackages
#

TRIBITS_PACKAGE_DECL( Teuchos ENABLE_SHADOWING_WARNINGS CLEANED )

SET(${PACKAGE_NAME}_VERSION "1.5d")

#
# B) Platform-specific checks
#

TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

INCLUDE(CheckForAbiDemangle)
CHECK_FOR_ABI_DEMANGLE(HAVE_GCC_ABI_DEMANGLE)

ASSERT_DEFINED(TPL_ENABLE_BinUtils CMAKE_CXX_COMPILER_ID)
IF(TPL_ENABLE_BinUtils)

  IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    SET(${PACKAGE_NAME}_ENABLE_STACKTRACE_DEFAULT ON)
  ELSE()
    SET(${PACKAGE_NAME}_ENABLE_STACKTRACE_DEFAULT OFF)
  ENDIF()

  TRIBITS_ADD_OPTION_AND_DEFINE(
    ${PACKAGE_NAME}_ENABLE_STACKTRACE
    HAVE_TEUCHOS_STACKTRACE
    "Enable Teuchos to show stack track for GCC (or Intel) on segfaults and uncaught exceptions."
    ${${PACKAGE_NAME}_ENABLE_STACKTRACE_DEFAULT}
    )
  
  IF (${PACKAGE_NAME}_ENABLE_STACKTRACE)
    SET(HAVE_TEUCHOS_LINK ON)
    SET(HAVE_TEUCHOS_BFD ON)
  ENDIF()

ELSE()

  SET(${PACKAGE_NAME}_ENABLE_STACKTRACE OFF)

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_STACKTRACE)
  TRIBITS_ADD_OPTION_AND_DEFINE(
    ${PACKAGE_NAME}_ENABLE_DEFAULT_STACKTRACE
    HAVE_TEUCHOS_DEFAULT_STACKTRACE
    "Enable Teuchos expansion of the stacktrace by default when exceptions are caught."
    ${${PACKAGE_NAME}_ENABLE_DEBUG}
    )
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_DEFAULT_STACKTRACE OFF)
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}Core_ENABLE_Boost)
IF (${PACKAGE_NAME}Core_ENABLE_Boost)
  INCLUDE(CheckForBoostIsPolymorphic)
  CHECK_FOR_BOOST_IS_POLYMORPHIC(HAVE_TEUCHOSCORE_BOOST_IS_POLYMORPHIC)
ENDIF()

ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}Numerics)
IF (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}Numerics)

  INCLUDE(CheckBlasFloat)
  CHECK_BLAS_FLOAT(HAVE_TEUCHOS_BLASFLOAT)
  IF(NOT HAVE_TEUCHOS_BLASFLOAT)
    # try apple veclib first
    IF(CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
      INCLUDE(CheckBlasFloatAppleVeclibBugfix)
      CHECK_BLAS_FLOAT_APPLE_VECLIB_BUGFIX(HAVE_TEUCHOS_BLASFLOAT_APPLE_VECLIB_BUGFIX)
      IF(HAVE_TEUCHOS_BLASFLOAT_APPLE_VECLIB_BUGFIX)
        GLOBAL_SET(HAVE_TEUCHOS_BLASFLOAT 1)
      ENDIF()
    ENDIF()
    # if veclib failed, general purpose solution
    IF(NOT HAVE_TEUCHOS_BLASFLOAT)
      INCLUDE(CheckBlasFloatReturnsDouble)
      CHECK_BLAS_FLOAT_DOUBLE_RETURN(HAVE_TEUCHOS_BLASFLOAT_DOUBLE_RETURN)
      IF(HAVE_TEUCHOS_BLASFLOAT_DOUBLE_RETURN)
        GLOBAL_SET(HAVE_TEUCHOS_BLASFLOAT 1)
      ENDIF()
    ENDIF()
  ENDIF()

  IF (HAVE_TEUCHOS_BLASFLOAT)
    INCLUDE(CheckSlapy2Problem)
    CHECK_SLAPY2_PROBLEM(HAVE_SLAPY2_PROBLEM)
    IF (HAVE_SLAPY2_PROBLEM) 
      INCLUDE(CheckSlapy2ReturnsDouble)
      CHECK_SLAPY2_RETURNS_DOUBLE(HAVE_SLAPY2_DOUBLE_RETURN)
      IF(NOT HAVE_SLAPY2_DOUBLE_RETURN)
        GLOBAL_SET(HAVE_TEUCHOS_BLASFLOAT OFF)
        MESSAGE(STATUS "SLAPY2 failure: setting HAVE_TEUCHOS_BLASFLOAT to OFF")
      ENDIF()
    ENDIF()
  ENDIF()

  INCLUDE(CheckLapackLARND)
  CHECK_LAPACK_LARND(HAVE_TEUCHOS_LAPACKLARND)

ENDIF()

SET( ${PACKAGE_NAME}_ORDINAL_TYPE ptrdiff_t
  CACHE STRING
  "Select the data type for Teuchos_Ordinal (Teuchos::Ordinal)." )

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

TRIBITS_SET_ST_FOR_DEV_MODE(ST_FOR_DEV_PT_FOR_RELEASE)

ASSERT_DEFINED(Trilinos_ENABLE_COMPLEX_DOUBLE)
ASSERT_DEFINED(Trilinos_ENABLE_COMPLEX_FLOAT)
IF (Trilinos_ENABLE_COMPLEX_DOUBLE OR Trilinos_ENABLE_COMPLEX_FLOAT)
  SET(${PACKAGE_NAME}_ENABLE_COMPLEX_DEFAULT  TRUE)
ELSE()
  # If Teuchos is configured as its own CMake project, independent of
  # Trilinos, then Teuchos_ENABLE_COMPLEX could be OFF by default.
  SET(${PACKAGE_NAME}_ENABLE_COMPLEX_DEFAULT  FALSE)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_COMPLEX
  HAVE_TEUCHOS_COMPLEX
  "Enable Teuchos features for complex data type."
  ${${PACKAGE_NAME}_ENABLE_COMPLEX_DEFAULT}
  )

# ToDo: Above, split ${PACKAGE_NAME}_ENABLE_COMPLEX into
# ${PACKAGE_NAME}_ENABLE_COMPLEX_FLOAT and
# ${PACKAGE_NAME}_ENABLE_COMPLEX_DOUBLE.  See Trilinos GitHub #362.

ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_CXX11)
SET(${PACKAGE_NAME}_ENABLE_LONG_LONG_INT_DEFAULT ${${PROJECT_NAME}_ENABLE_CXX11})
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_LONG_LONG_INT
  HAVE_TEUCHOS_LONG_LONG_INT
  "Enable Teuchos features for the 'long long' data type.  This is ON by default if and only if C++11 is enabled.  Please remember that 'long long' only entered the C++ language standard with C++11, even though many compilers supported it before."
  ${${PACKAGE_NAME}_ENABLE_LONG_LONG_INT_DEFAULT}
  )

IF(WIN32 AND NOT CYGWIN AND NOT ${PACKAGE_NAME}_ENABLE_LONG_LONG_INT)
  # Can only provide __int64 support if it will not clash with 
  # long long int (see Trilinos bug 5801).
  INCLUDE(CheckFor__int64)
  CHECK_FOR___INT64(HAVE_TEUCHOS___INT64)
ENDIF()

# Check whether the C++ compiler can compile __attribute__((constructor)).
INCLUDE(CheckCXXAttributeConstructorCompiles)
CHECK_CXX_ATTRIBUTE_CONSTRUCTOR_COMPILES(HAVE_TEUCHOS_CXX_ATTRIBUTE_CONSTRUCTOR)
# Check whether the C++ compiler can compile __attribute__((weak)).
INCLUDE(CheckCXXAttributeWeakCompiles)
CHECK_CXX_ATTRIBUTE_WEAK_COMPILES(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)
# Some compilers don't implement __attribute__((weak)), but they might
# implement "#pragma weak" instead.  Check this.
INCLUDE(CheckCXXPragmaWeakCompiles)
CHECK_CXX_PRAGMA_WEAK_COMPILES(HAVE_TEUCHOS_CXX_PRAGMA_WEAK)
# Check whether building dynamic shared libraries or static libraries.
# If building dynamic shared libraries, define a macro.  In that case,
# fixing Bug 6392 doesn't need any of the three things above.
GLOBAL_SET(HAVE_TEUCHOS_DYNAMIC_LIBS OFF) # to revise below
IF (DEFINED BUILD_SHARED_LIBS)
  # Don't check the value unless we know that it's defined.
  IF (BUILD_SHARED_LIBS)
    GLOBAL_SET(HAVE_TEUCHOS_DYNAMIC_LIBS ON)
  ENDIF ()
ENDIF ()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_THREAD_SAFE
  HAVE_TEUCHOS_THREAD_SAFE
  "Enable thread safe code for RCP classes."
  ${Trilinos_ENABLE_THREAD_SAFE}
  )

IF (${PACKAGE_NAME}_ENABLE_THREAD_SAFE AND NOT ${PROJECT_NAME}_ENABLE_CXX11)
  MESSAGE(FATAL_ERROR
    "You set ${PACKAGE_NAME}_ENABLE_THREAD_SAFE=ON, but ${PACKAGE_NAME}' support"
    " for CXX11 is not enabled (${PROJECT_NAME}_ENABLE_CXX11=OFF)."
    "  This is not allowed.  Please enable ${PACKAGE_NAME}_ENABLE_CXX11 in"
    " ${PACKAGE_NAME} before attempting to enable ${PACKAGE_NAME}_ENABLE_THREAD_SAFE"
    " or leave ${PACKAGE_NAME}_ENABLE_THREAD_SAFE off.")
ENDIF ()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_THREAD_SAFE)
ASSERT_DEFINED(${PACKAGE_NAME}Core_ENABLE_Pthread)
COMBINED_OPTION(
  ${PACKAGE_NAME}_ENABLE_THREAD_SAFE_UNIT_TESTS
  DEP_OPTIONS_NAMES
    ${PACKAGE_NAME}_ENABLE_THREAD_SAFE
    ${PACKAGE_NAME}Core_ENABLE_Pthread
  DOCSTR
    "Enable thread safe unit tests for ${PACKAGE_NAME} Memory Management classes."
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_ABC
  HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  "Enable ${PACKAGE_NAME} features for array bounds checking."
  ${${PROJECT_NAME}_ENABLE_DEBUG}
  )

IF (${PACKAGE_NAME}_ENABLE_COMPLEX)

  INCLUDE(CheckCXXComplexBlasProblem)
  CHECK_CXX_COMPLEX_BLAS_PROBLEM(HAVE_COMPLEX_BLAS_PROBLEM)

  IF (HAVE_COMPLEX_BLAS_PROBLEM)
    INCLUDE(CheckCXXComplexBlasProblemCanBeFixed)
    IF(HAVE_TEUCHOS_BLASFLOAT_APPLE_VECLIB_BUGFIX)
      CHECK_COMPLEX_BLAS_VECLIB_OKAY(HAVE_VECLIB_COMPLEX_BLAS)
    ELSE()
      CHECK_CXX_COMPLEX_BLAS_PROBLEM_CAN_BE_FIXED(HAVE_FIXABLE_COMPLEX_BLAS_PROBLEM)
    ENDIF()
  ENDIF()
  GLOBAL_SET(HAVE_COMPLEX_BLAS 
             NOT HAVE_COMPLEX_BLAS_PROBLEM OR HAVE_FIXABLE_COMPLEX_BLAS_PROBLEM OR HAVE_VECLIB_COMPLEX_BLAS)

ELSE() # NOT ${PACKAGE_NAME}_ENABLE_COMPLEX

  # Forbid users from setting HAVE_COMPLEX_BLAS=ON, if
  # Teuchos_ENABLE_COMPLEX is OFF.  Use two separate IF statements,
  # because "DEFINED X AND NOT X" could be bad syntax if X is NOT
  # defined.
  IF (DEFINED HAVE_COMPLEX_BLAS)
    IF (NOT HAVE_COMPLEX_BLAS)
      MESSAGE(FATAL_ERROR "Teuchos: ${PACKAGE_NAME}_ENABLE_COMPLEX is OFF, yet you set HAVE_COMPLEX_BLAS=ON.  Please do not do this.  You normally do not need to set HAVE_COMPLEX_BLAS explicitly.  If for some reason you do, then it must have the same value as ${PACKAGE_NAME}_ENABLE_COMPLEX.")
    ENDIF()
  ENDIF()

ENDIF() # ${PACKAGE_NAME}_ENABLE_COMPLEX

#
# C) Set up package-specific options
#

IF (Trilinos_ENABLE_FLOAT)
  SET(${PACKAGE_NAME}_ENABLE_FLOAT_DEFAULT  TRUE)
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_FLOAT_DEFAULT  FALSE)
  # ABove will be set in case Teuchos is confiugred as its own CMake project!
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_FLOAT
  HAVE_TEUCHOS_FLOAT
  "Enable Teuchos features for float data type."
  ${${PACKAGE_NAME}_ENABLE_FLOAT_DEFAULT}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_DEBUG_RCP_NODE_TRACING
  HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING
  "At compile time, enable debug-mode RCP node tracing which can't be turned off."
  ${${PACKAGE_NAME}_ENABLE_DEBUG}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_EXTENDED
  HAVE_TEUCHOS_EXTENDED
  "Enable extended Teuchos features."
  ON
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_C_EXCEPTIONS
  HAVE_TEUCHOS_C_EXCEPTIONS
  "Enable support for the C macro TEUCHOS_EXIT(...) to thrown an exception."
  OFF
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_GCC_DEMANGLE
  HAVE_TEUCHOS_DEMANGLE
  "Enable Teuchos calls to GCC demangle library function if it exists (see HAVE_TEUCHOS_BLASFLOAT)."
  ON
  )

IF (Teuchos_ENABLE_QD)
  TRIBITS_ADD_OPTION_AND_DEFINE(
    ${PACKAGE_NAME}_PROMOTE_DOUBLE_VIA_QD
    HAVE_TEUCHOS_DOUBLE_TO_QD
    "Instruct Teuchos::ScalarTraits to promote double precision to the double-double of the QD library (see HAVE_TEUCHOS_QD)."
    NO
    )
ENDIF()

IF (Teuchos_ENABLE_ARPREC)
  TRIBITS_ADD_OPTION_AND_DEFINE(
    ${PACKAGE_NAME}_PROMOTE_DOUBLE_VIA_ARPREC
    HAVE_TEUCHOS_DOUBLE_TO_ARPREC
    "Instruct Teuchos::ScalarTraits to promote double precision to the double-double of the QD library (see HAVE_TEUCHOS_ARPREC)."
    NO
    )
ENDIF()

IF (HAVE_TEUCHOS_DOUBLE_TO_ARPREC AND HAVE_TEUCHOS_DOUBLE_TO_QD) 
    MESSAGE("Cannot promote double to both ARPREC and QD. Will use ARPREC.")
    SET(HAVE_TEUCHOS_DOUBLE_TO_QD OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_GLOBALLY_REDUCE_UNITTEST_RESULTS
  HAVE_TEUCHOS_GLOBALLY_REDUCE_UNITTEST_RESULTS
  "Set the Teuchos::UnitTestRepository to globally reduce test pass/fail results across all MPI ranks."
  OFF
  )

ASSERT_DEFINED(TPL_ENABLE_Valgrind)
IF(TPL_ENABLE_Valgrind)
  # Heap profiling using massif
  TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_TIME_MASSIF_SNAPSHOTS
    HAVE_TEUCHOS_TIME_MASSIF_SNAPSHOTS
    "Enable memory profiling using massif."
    OFF)
  IF (${PACKAGE_NAME}_TIME_MASSIF_SNAPSHOTS)
    MESSAGE(STATUS "Enabling memory profiling using massif")
  ENDIF()
ENDIF()

# Enabling Kokkos profiling hooks in Teuchos timers requires that the
# Kokkos package be enabled, and that KOKKOS_ENABLE_PROFILING is ON.

IF(DEFINED ${PROJECT_NAME}_ENABLE_Kokkos AND ${PROJECT_NAME}_ENABLE_Kokkos)
  ASSERT_DEFINED(KOKKOS_ENABLE_PROFILING)
  SET(${PACKAGE_NAME}_KOKKOS_PROFILING_DEFAULT ${KOKKOS_ENABLE_PROFILING})
ELSE()
  SET(${PACKAGE_NAME}_KOKKOS_PROFILING_DEFAULT OFF)
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_KOKKOS_PROFILING_DEFAULT)
TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_KOKKOS_PROFILING
  HAVE_TEUCHOS_KOKKOS_PROFILING
  "Wrap every Teuchos timer with a Kokkos profiling region."
  "${${PACKAGE_NAME}_KOKKOS_PROFILING_DEFAULT}")
IF (${PACKAGE_NAME}_KOKKOS_PROFILING)
  MESSAGE(STATUS "Wrapping every Teuchos timer with a Kokkos profiling region.")
ENDIF()


#
# D) Process the subpackages for Teuchos
#

# For use by other packages as well as Teuchos tests
GLOBAL_SET(TEUCHOS_STD_UNIT_TEST_MAIN
  "${PACKAGE_SOURCE_DIR}/core/test/UnitTest/Teuchos_StandardUnitTestMain.cpp")

# For use by other packages as well as Teuchos tests
IF (${PROJECT_NAME}_ENABLE_TeuchosComm)
  GLOBAL_SET(TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN
    "${PACKAGE_SOURCE_DIR}/comm/test/UnitTesting/Teuchos_StandardParallelUnitTestMain.cpp")
ENDIF()

TRIBITS_PROCESS_SUBPACKAGES()

TRIBITS_PACKAGE_DEF()

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_EXCLUDE_FILES(
  config.h.in
  doc/images
  example/config.h.in
  )

#
# E) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()


ELSE(TRIBITS_PROCESSING_PACKAGE)


#
# Processing as a TriBITS project
#

# This CMakeLists.txt file is being processed as the TriBITS projects's base
# CMakeLists.txt file!  (See comments at bottom of this file.)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11 FATAL_ERROR)
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/ProjectName.cmake")
PROJECT(${PROJECT_NAME} NONE)
SET(${PROJECT_NAME}_TRIBITS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/tribits"
  CACHE PATH "By default assume Teuchos is in Trilinos")
INCLUDE("${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake")
SET(TEUCHOS_STANDALONE_PACKAGE TRUE)
SET(${PROJECT_NAME}_ENABLE_CXX11_DEFAULT ON)
SET(${PROJECT_NAME}_ENABLE_Teuchos  ON  CACHE  BOOL
  "Set by default in teuchos/CMakeLists.txt")
TRIBITS_PROJECT()


ENDIF(TRIBITS_PROCESSING_PACKAGE)
