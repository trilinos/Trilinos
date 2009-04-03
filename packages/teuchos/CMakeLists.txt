INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

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

INCLUDE(CheckBlasFloat)
CHECK_BLAS_FLOAT(HAVE_TEUCHOS_BLASFLOAT)

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_COMPLEX
  HAVE_TEUCHOS_COMPLEX
  "Enable Teuchos features for complex data type."
  ON
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_LONG_LONG
  HAVE_TEUCHOS_LONG_LONG_INT
  "Enable Teuchos features for long long data type."
  OFF
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
  ON
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_DEBUG
  HAVE_TEUCHOS_DEBUG
  "Enable a host of runtime debug checking."
  ${Trilinos_ENABLE_DEBUG}
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION
  HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
  "Turn on explicit instantiation in Teuchos."
  ${Trilinos_ENABLE_EXPLICIT_INSTANTIATION}
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
