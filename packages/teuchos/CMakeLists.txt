INCLUDE(PackageMacros)

#
# A) Define the package
#

PACKAGE( Teuchos ENABLE_SHADOWING_WARNINGS CLEANED )

SET(${PACKAGE_NAME}_VERSION "1.5d")

#
# B) Platform-specific checks
#

INCLUDE(CheckForAbiDemangle)
CHECK_FOR_ABI_DEMANGLE(HAVE_GCC_ABI_DEMANGLE)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Boost)
IF (${PACKAGE_NAME}_ENABLE_Boost)
  INCLUDE(CheckForBoostIsPolymorphic)
  CHECK_FOR_BOOST_IS_POLYMORPHIC(HAS_TEUCHOS_BOOST_IS_POLYMORPHIC)
ENDIF()

INCLUDE(CheckBlasFloat)
CHECK_BLAS_FLOAT(HAVE_TEUCHOS_BLASFLOAT)

SET( ${PACKAGE_NAME}_ORDINAL_TYPE ptrdiff_t
  CACHE STRING
  "Select the data type for Teuchos_Ordinal (Teuchos::Ordinal)." )

PACKAGE_ADD_DEBUG_OPTION()

PACKAGE_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

PACKAGE_ADD_EXPLICIT_INSTANTIATION_OPTION()

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_COMPLEX
  HAVE_TEUCHOS_COMPLEX
  "Enable Teuchos features for complex data type."
  ON
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_LONG_LONG_INT
  HAVE_TEUCHOS_LONG_LONG_INT
  "Enable Teuchos features for long long data type."
  OFF
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_ABC
  HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  "Enable Teuchos features for array bounds checking."
  ${Trilinos_ENABLE_DEBUG}
  )

IF (${PACKAGE_NAME}_ENABLE_COMPLEX)

  INCLUDE(CheckCXXComplexBlasProblem)
  CHECK_CXX_COMPLEX_BLAS_PROBLEM(HAVE_COMPLEX_BLAS_PROBLEM)
  
  IF (HAVE_COMPLEX_BLAS_PROBLEM)
    INCLUDE(CheckCXXComplexBlasProblemCanBeFixed)
    CHECK_CXX_COMPLEX_BLAS_PROBLEM_CAN_BE_FIXED(
      HAVE_FIXABLE_COMPLEX_BLAS_PROBLEM)
  ENDIF()

ENDIF()

#
# C) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_FLOAT
  HAVE_TEUCHOS_FLOAT
  "Enable Teuchos features for float data type."
  ${${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE}
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_DEBUG_RCP_NODE_TRACING
  HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING
  "At compile time, enable debug-mode RCP node tracing which can't be turned off."
  ${${PACKAGE_NAME}_ENABLE_DEBUG}
  )

OPTION(
  ${PACKAGE_NAME}_ENABLE_PERFORMANCE_TESTS
  "Override enable ${PACKAGE_NAME} performance tests."
  ${Trilinos_ENABLE_PERFORMANCE_TESTS}
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_EXTENDED
  HAVE_TEUCHOS_EXTENDED
  "Enable extended Teuchos features."
  ON
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_C_EXCEPTIONS
  HAVE_TEUCHOS_C_EXCEPTIONS
  "Enable support for the C macro TEUCHOS_EXIT(...) to thrown an exception."
  OFF
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_GCC_DEMANGLE
  HAVE_TEUCHOS_DEMANGLE
  "Enable Teuchos calls to GCC demangle library function if it exists (see HAVE_TEUCHOS_BLASFLOAT)."
  ON
  )

IF (Teuchos_ENABLE_QD)
  ADD_OPTION_AND_DEFINE(
    ${PACKAGE_NAME}_PROMOTE_DOUBLE_VIA_QD
    HAVE_TEUCHOS_DOUBLE_TO_QD
    "Instruct Teuchos::ScalarTraits to promote double precision to the double-double of the QD library (see HAVE_TEUCHOS_QD)."
    NO
    )
ENDIF()

IF (Teuchos_ENABLE_ARPREC)
  ADD_OPTION_AND_DEFINE(
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

#
# D) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

PACKAGE_ADD_PERFORMANCE_TEST_DIRECTORIES(performance)

#
# E) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
