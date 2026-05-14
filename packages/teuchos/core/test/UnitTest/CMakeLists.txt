#
# UnitTest_UnitTests
#


TRIBITS_ADD_EXECUTABLE(
  UnitTest_UnitTests
  SOURCES
    UnitTest_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  )


TRIBITS_ADD_TEST(
  UnitTest_UnitTests
  POSTFIX_AND_ARGS_0 test_names --details=TEST_NAMES
  POSTFIX_AND_ARGS_1 final_results --details=FINAL_RESULTS
  POSTFIX_AND_ARGS_2 all --details=ALL
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


#
# Simple_UnitTests
#


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Simple_UnitTests
  SOURCES
    Int_UnitTests.cpp
    PrintDouble_UnitTests.cpp
    TemplateFunc_UnitTests.cpp
    vector_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

# PASS_REGULAR_EXPRESSION below must mirror the binary's "Summary"
# line byte-for-byte. The binary's test count is fixed at compile time by
# HAVE_TEUCHOS_INST_* (see Teuchos_UnitTestHelpers.hpp), not by the higher-
# level ${PARENT_PACKAGE_NAME}_ENABLE_{COMPLEX,FLOAT}, from which the INST_*
# values are only derived as defaults. The two diverge whenever the user
# overrides Trilinos_ENABLE_COMPLEX_{FLOAT,DOUBLE} or ${PARENT_PACKAGE_NAME}_INST_*
# directly (see packages/teuchos/CMakeLists.txt). Concrete trap: the default
# Teuchos_INST_COMPLEX_FLOAT for (ENABLE_COMPLEX=ON, ENABLE_FLOAT=OFF) is OFF,
# so an _ENABLE_*-based prediction overcounts by one for that very common
# configuration. Predicting from _INST_* makes the divergence impossible by
# construction.
#
# Base count = 16, justified by the unconditional registrations:
#   Int_UnitTests.cpp           2  Basic, Assignment
#   PrintDouble_UnitTests.cpp   1  Basic
#   vector_UnitTests.cpp       12  4 tests * {int, float, double}; the type
#                                  set is hard-coded and ignores *_INST_FLOAT
#   TemplateFunc_UnitTests.cpp  1  double leg of
#                                  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES,
#                                  which has no HAVE_TEUCHOS_INST_DOUBLE guard
set(Simple_UnitTests_total 16)
if (${PARENT_PACKAGE_NAME}_INST_FLOAT)
  math(EXPR Simple_UnitTests_total "${Simple_UnitTests_total} + 1")
endif()
if (${PARENT_PACKAGE_NAME}_INST_COMPLEX_FLOAT)
  math(EXPR Simple_UnitTests_total "${Simple_UnitTests_total} + 1")
endif()
if (${PARENT_PACKAGE_NAME}_INST_COMPLEX_DOUBLE)
  math(EXPR Simple_UnitTests_total "${Simple_UnitTests_total} + 1")
endif()

# run = total - 1: the minus_one variant filters exactly vector<float>::constAt
# (see ARGS below), which is part of the unconditional base above and is
# therefore guaranteed to be present in every supported configuration.
math(EXPR Simple_UnitTests_run "${Simple_UnitTests_total} - 1")
set(PASS_STRING
  "Summary: total = ${Simple_UnitTests_total}, run = ${Simple_UnitTests_run}, passed = ${Simple_UnitTests_run}, failed = 0")


TRIBITS_ADD_TEST(
  Simple_UnitTests
  NAME_POSTFIX minus_one
  ARGS "--not-unit-test=vector_float_constAt_UnitTest"
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "${PASS_STRING}"
  )


#
# UnitTest_BadUnitTest
#


TRIBITS_ADD_EXECUTABLE(
  UnitTest_BadUnitTest
  COMM serial mpi  
  SOURCES
    Int_UnitTests.cpp
    TemplateFunc_UnitTests.cpp
    Failing_UnitTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  )


TRIBITS_ADD_TEST(
  UnitTest_BadUnitTest
  NAME UnitTest_BadUnitTest_final_results
  ARGS "--show-test-details=FINAL_RESULTS"
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "Int_BadAssignment_UnitTest"
  )


TRIBITS_ADD_TEST(
  UnitTest_BadUnitTest
  NAME UnitTest_BadUnitTest_end_result_failed
  ARGS "--teuchos-suppress-startup-banner"
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "End Result: TEST FAILED"
  )


TRIBITS_ADD_TEST(
  UnitTest_BadUnitTest
  NAME UnitTest_BadUnitTest_end_result_totals
  ARGS "--teuchos-suppress-startup-banner"
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION ", failed = 2"
  )


TRIBITS_ADD_TEST(
  UnitTest_BadUnitTest
  NAME UnitTest_BadUnitTest_Int_BadAssignment_failed
  ARGS
    "--show-test-details=FINAL_RESULTS --group-name=Int --test-name=BadAssignment"
    "--show-test-details=TEST_NAMES --group-name=Int --test-name=BadAssignment"
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "[FAILED]"
  )


#
# Hashtable unit tests
#


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Hashtable_UnitTest
  SOURCES
    Hashtable_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


#
# Range1D unit tests
#


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Range1D_UnitTests
  SOURCES
    Range1D_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


#
# VerbosityLevel unit tests
#

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  VerbosityLevel_UnitTests
  SOURCES
    VerbosityLevel_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


#
# GlobalMPISession unit tests
#


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  GlobalMPISessionUnitTests
  SOURCES GlobalMPISession_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
TRIBITS_ADD_TEST(
  GlobalMPISessionUnitTests
  COMM mpi
  STANDARD_PASS_OUTPUT
  )


#
# SystemInformation unit tests
#

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SystemInformation_UnitTests
  SOURCES
    SystemInformation_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
