# -*- cmake -*-

INCLUDE(TribitsAddExecutableAndTest)
#INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Unit test the Slice class
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Slice_UnitTests
  SOURCES
    Slice_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS --teuchos-suppress-startup-banner
  STANDARD_PASS_OUTPUT
  )

# Unit test the MDArray class
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MDArray_UnitTests
  SOURCES
    MDArray_UnitTests.cpp
    MDArray_UnitTest_helpers.hpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS --teuchos-suppress-startup-banner
  STANDARD_PASS_OUTPUT
  )

# Unit test the MDArrayView class
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MDArrayView_UnitTests
  SOURCES
    MDArrayView_UnitTests.cpp
    MDArray_UnitTest_helpers.hpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS --teuchos-suppress-startup-banner
  STANDARD_PASS_OUTPUT
  )

# Unit test the MDArrayRCP class
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MDArrayRCP_UnitTests
  SOURCES
    MDArrayRCP_UnitTests.cpp
    MDArray_UnitTest_helpers.hpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS --teuchos-suppress-startup-banner
  STANDARD_PASS_OUTPUT
  )

# Performance test the MDArray class
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MDArray_PerformanceTests
  NAME_POSTFIX basic
  CATEGORIES BASIC PERFORMANCE
  SOURCES
    MDArray_Performance_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

#TRIBITS_COPY_FILES_TO_BINARY_DIR(target_name
#  SOURCE_FILES source_files
#  EXEDEPS exedeps
#  )
