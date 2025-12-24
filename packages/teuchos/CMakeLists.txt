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
# A) Forward declare the package so that certain options are also defined for
# subpackages
#

TRIBITS_PACKAGE_DECL( Teuchos ENABLE_SHADOWING_WARNINGS CLEANED )

SET(${PACKAGE_NAME}_VERSION "1.5d")

#
# B) Platform-specific checks
#

TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

#
# The check needs to run code, so give cross-compiling users a bypass.
#
IF (DEFINED HAVE_GCC_ABI_DEMANGLE)
  MESSAGE(STATUS "You set HAVE_GCC_ABI_DEMANGLE=${HAVE_GCC_ABI_DEMANGLE}; this suggests you are cross compiling.  If you're not cross compiling, best practice is not to set this; Trilinos will figure it out for you.")
ELSE()
  INCLUDE(CheckForAbiDemangle)
  CHECK_FOR_ABI_DEMANGLE(HAVE_GCC_ABI_DEMANGLE)
ENDIF()
ASSERT_DEFINED(HAVE_GCC_ABI_DEMANGLE)

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
  #
  # The check needs to run code, so give cross-compiling users a bypass.
  #
  IF (DEFINED HAVE_TEUCHOS_BLASFLOAT)
    MESSAGE(STATUS "You set HAVE_TEUCHOS_BLASFLOAT=${HAVE_TEUCHOS_BLASFLOAT}; this suggests you are cross compiling.  If you're not cross compiling, best practice is not to set this; Trilinos will figure it out for you.")
  ELSE()
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
  ENDIF()
  ASSERT_DEFINED(HAVE_TEUCHOS_BLASFLOAT)

  #
  # The check needs to run code, so give cross-compiling users a bypass.
  #
  IF (DEFINED HAVE_TEUCHOS_LAPACKLARND)
    MESSAGE(STATUS "You set HAVE_TEUCHOS_LAPACKLARND=${HAVE_TEUCHOS_LAPACKLARND}; this suggests you are cross compiling.  If you're not cross compiling, best practice is not to set this; Trilinos will figure it out for you.")
  ELSE()
    INCLUDE(CheckLapackLARND)
    CHECK_LAPACK_LARND(HAVE_TEUCHOS_LAPACKLARND)
  ENDIF()
  ASSERT_DEFINED(HAVE_TEUCHOS_LAPACKLARND)
ENDIF()

SET( ${PACKAGE_NAME}_ORDINAL_TYPE ptrdiff_t
  CACHE STRING
  "Select the data type for Teuchos_Ordinal (Teuchos::Ordinal)." )

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

TRIBITS_SET_ST_FOR_DEV_MODE(ST_FOR_DEV_PT_FOR_RELEASE)

#
# Set cache vars for support/instantiations of various types
#
# NOTE: The else() statements for the default are for the case where Teuchos
# is configured as its own CMake project and variables tarting with Trilinos_
# have not been set.
#

IF (Trilinos_ENABLE_LONG_DOUBLE)
  SET(${PACKAGE_NAME}_ENABLE_LONG_DOUBLE_DEFAULT  TRUE)
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_LONG_DOUBLE_DEFAULT  FALSE)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_LONG_DOUBLE
  HAVE_TEUCHOS_LONG_DOUBLE
  "Enable Teuchos features for long double data type."
  ${${PACKAGE_NAME}_ENABLE_LONG_DOUBLE_DEFAULT}
  )

IF (Trilinos_ENABLE_FLOAT)
  SET(${PACKAGE_NAME}_ENABLE_FLOAT_DEFAULT  TRUE)
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_FLOAT_DEFAULT  FALSE)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_FLOAT
  HAVE_TEUCHOS_FLOAT
  "Enable Teuchos features for float data type."
  ${${PACKAGE_NAME}_ENABLE_FLOAT_DEFAULT}
  )

IF (Trilinos_ENABLE_COMPLEX_DOUBLE OR Trilinos_ENABLE_COMPLEX_FLOAT)
  SET(${PACKAGE_NAME}_ENABLE_COMPLEX_DEFAULT  TRUE)
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_COMPLEX_DEFAULT  FALSE)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_COMPLEX
  HAVE_TEUCHOS_COMPLEX
  "Enable basic Teuchos features for complex data types."
  ${${PACKAGE_NAME}_ENABLE_COMPLEX_DEFAULT}
  )

IF (${PACKAGE_NAME}_ENABLE_FLOAT)
  SET(${PACKAGE_NAME}_INST_FLOAT_DEFAULT  TRUE)
ELSE()
  SET(${PACKAGE_NAME}_INST_FLOAT_DEFAULT  FALSE)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_FLOAT
  HAVE_TEUCHOS_INST_FLOAT
  "Do explicit instantaitions for float scalar type."
  ${${PACKAGE_NAME}_INST_FLOAT_DEFAULT}
  )

IF (NOT Trilinos_ENABLE_COMPLEX_FLOAT STREQUAL "")
  SET(${PACKAGE_NAME}_INST_COMPLEX_FLOAT_DEFAULT ${Trilinos_ENABLE_COMPLEX_FLOAT})
ELSEIF (${PACKAGE_NAME}_ENABLE_COMPLEX AND ${PACKAGE_NAME}_ENABLE_FLOAT)
  SET(${PACKAGE_NAME}_INST_COMPLEX_FLOAT_DEFAULT  TRUE)
ELSE()
  SET(${PACKAGE_NAME}_INST_COMPLEX_FLOAT_DEFAULT  FALSE)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_COMPLEX_FLOAT
  HAVE_TEUCHOS_INST_COMPLEX_FLOAT
  "Do explicit instantiations for complex<float> scalar type."
  ${${PACKAGE_NAME}_INST_COMPLEX_FLOAT_DEFAULT}
  )

IF (Trilinos_ENABLE_COMPLEX_DOUBLE OR ${PACKAGE_NAME}_ENABLE_COMPLEX)
  SET(${PACKAGE_NAME}_INST_COMPLEX_DOUBLE_DEFAULT  TRUE)
ELSE()
  SET(${PACKAGE_NAME}_INST_COMPLEX_DOUBLE_DEFAULT  FALSE)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_COMPLEX_DOUBLE
  HAVE_TEUCHOS_INST_COMPLEX_DOUBLE
  "Do explicit instantiations for complex<float> scalar type."
  ${${PACKAGE_NAME}_INST_COMPLEX_DOUBLE_DEFAULT}
  )

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


# As per the discussion in #11921, we explicitly enable thread safe RCPs,
# if Kokkos has a host thread-parallel backend
SET(Teuchos_ENABLE_THREAD_SAFE_Default ${Trilinos_ENABLE_THREAD_SAFE})
IF(DEFINED ${PROJECT_NAME}_ENABLE_Kokkos AND ${PROJECT_NAME}_ENABLE_Kokkos)
  IF( (DEFINED Kokkos_ENABLE_THREADS AND Kokkos_ENABLE_THREADS) OR
      (DEFINED Kokkos_ENABLE_OPENMP AND Kokkos_ENABLE_OPENMP) )
    IF (NOT ${PACKAGE_NAME}_ENABLE_THREAD_SAFE)
      MESSAGE(WARNING "A Kokkos host thread-parallel backend is enabled, requiring Trilinos to enable thread-safe RCP classes.  You can do this manually by specifying Teuchos_ENABLE_THREAD_SAFE=ON")
      SET(Teuchos_ENABLE_THREAD_SAFE_Default ON)
    ENDIF()
  ENDIF()
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_THREAD_SAFE
  HAVE_TEUCHOS_THREAD_SAFE
  "Enable thread safe code for RCP classes."
  ${Teuchos_ENABLE_THREAD_SAFE_Default}
  )



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
  #
  # The check needs to run code, so give cross-compiling users a
  # bypass.  This is not a complete check, but we have an error test
  # below that explains the required settings.  (In the cross
  # compilation case, running code here could get the wrong answer.
  # Thus, it's better to error out than to continue.)
  #
  IF (DEFINED HAVE_COMPLEX_BLAS_PROBLEM)
    MESSAGE(STATUS "You set HAVE_COMPLEX_BLAS_PROBLEM=${HAVE_COMPLEX_BLAS_PROBLEM}; this suggests you are cross compiling.  If you're not cross compiling, best practice is not to set this; Trilinos will figure it out for you.")
  ELSE()
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
  ENDIF()
  ASSERT_DEFINED(HAVE_COMPLEX_BLAS_PROBLEM)
  IF (HAVE_COMPLEX_BLAS_PROBLEM)
    IF ((NOT (DEFINED HAVE_VECLIB_COMPLEX_BLAS)) AND
        (NOT (DEFINED HAVE_FIXABLE_COMPLEX_BLAS_PROBLEM)))
      MESSAGE(FATAL_ERROR "IF HAVE_COMPLEX_BLAS_PROBLEM is defined and TRUE, then you must define at least one of HAVE_VECLIB_COMPLEX_BLAS or HAVE_FIXABLE_COMPLEX_BLAS_PROBLEM.  This is most likely to happen if you are cross compiling and setting these options yourself, rather than relying on CMake to set the options.")
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

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_MODIFY_DEFAULTS_DURING_VALIDATION
  HAVE_TEUCHOS_MODIFY_DEFAULTS_DURING_VALIDATION
  "Modify default parameters during validation."
  ON
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
  ${PACKAGE_NAME}_ENABLE_STACKED_TIMER_IN_TIME_MONITOR
  HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
  "Adds timings from TimeMonitor to the StackedTimer object."
  ON
  )

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

#Option Teuchos_KOKKOS_PROFILING
IF(DEFINED ${PROJECT_NAME}_ENABLE_Kokkos AND ${PROJECT_NAME}_ENABLE_Kokkos)
  SET(${PACKAGE_NAME}_KOKKOS_PROFILING_DEFAULT ON)
ELSE()
  SET(${PACKAGE_NAME}_KOKKOS_PROFILING_DEFAULT OFF)
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_KOKKOS_PROFILING_DEFAULT)
TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_KOKKOS_PROFILING
  HAVE_TEUCHOS_KOKKOS_PROFILING
  "Wrap every Teuchos timer with a Kokkos profiling region."
  "${${PACKAGE_NAME}_KOKKOS_PROFILING_DEFAULT}")
IF (${PACKAGE_NAME}_KOKKOS_PROFILING)
  # Enabling Kokkos profiling hooks in Teuchos timers requires that the
  # Kokkos package be enabled
  IF(NOT(DEFINED ${PROJECT_NAME}_ENABLE_Kokkos AND ${PROJECT_NAME}_ENABLE_Kokkos))
    MESSAGE(FATAL_ERROR "Setting Teuchos_KOKKOS_PROFILING=ON requires that the Kokkos package is enabled.")
  ENDIF()
  MESSAGE(STATUS "Wrapping every Teuchos timer with a Kokkos profiling region.")
ENDIF()

#Option Teuchos_TIMER_KOKKOS_FENCE
SET(Teuchos_TIMER_KOKKOS_FENCE_DEFAULT OFF)
TRIBITS_ADD_OPTION_AND_DEFINE(
  Teuchos_TIMER_KOKKOS_FENCE
  HAVE_TEUCHOS_TIMER_KOKKOS_FENCE
  "Call Kokkos::fence() whenever a Teuchos timer starts or stops."
  "${Teuchos_TIMER_KOKKOS_FENCE_DEFAULT}")
IF (Teuchos_TIMER_KOKKOS_FENCE)
  # Enabling Kokkos::fence() in Teuchos timers requires that the Kokkos package be enabled
  IF(NOT(DEFINED ${PROJECT_NAME}_ENABLE_Kokkos AND ${PROJECT_NAME}_ENABLE_Kokkos))
    MESSAGE(FATAL_ERROR "Setting Teuchos_TIMER_KOKKOS_FENCE=ON requires that the Kokkos package is enabled.")
  ENDIF()
  MESSAGE(STATUS "Will call Kokkos::fence() when a Teuchos timer starts or stops.")
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
CMAKE_MINIMUM_REQUIRED(VERSION 3.10.0 FATAL_ERROR)
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/ProjectName.cmake")
PROJECT(${PROJECT_NAME} NONE)
SET(${PROJECT_NAME}_TRIBITS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/tribits"
  CACHE PATH "By default assume Teuchos is in Trilinos")
INCLUDE("${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake")
SET(TEUCHOS_STANDALONE_PACKAGE TRUE)
SET(${PROJECT_NAME}_ENABLE_Teuchos  ON  CACHE  BOOL
  "Set by default in teuchos/CMakeLists.txt")
TRIBITS_PROJECT()


ENDIF(TRIBITS_PROCESSING_PACKAGE)
