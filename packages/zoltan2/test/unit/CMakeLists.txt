INCLUDE(TribitsAddExecutableAndTest)

# TODO
# We want to emit "NOTRUN" when a test is skipped
# for some reason. How do we tell cmake?


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS
  AdaptersForTests.hpp
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  IdentifierMap
  SOURCES IdentifierMap.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  IdentifierTraits
  SOURCES IdentifierTraits.cpp 
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  XpetraTraits
  SOURCES XpetraTraits.cpp UserInputForTests.hpp ErrorHandlingForTests.hpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  AlltoAll
  SOURCES AlltoAll.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CommTest
  SOURCES CommTest.cpp
  NUM_MPI_PROCS 4
  COMM mpi serial
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Parameters
  SOURCES Parameters.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )
 
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  XpetraMatrixAdapter
  SOURCES XpetraMatrixAdapter.cpp UserInputForTests.hpp ErrorHandlingForTests.hpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )
 
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  XpetraGraphAdapter
  SOURCES XpetraGraphAdapter.cpp UserInputForTests.hpp ErrorHandlingForTests.hpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  XpetraVectorAdapter
  SOURCES XpetraVectorAdapter.cpp UserInputForTests.hpp ErrorHandlingForTests.hpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  XpetraMultiVectorAdapter
  SOURCES XpetraMultiVectorAdapter.cpp UserInputForTests.hpp ErrorHandlingForTests.hpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  GraphModel
  SOURCES GraphModel.cpp AdaptersForTests.hpp ErrorHandlingForTests.hpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  IdentifierModel
  SOURCES IdentifierModel.cpp AdaptersForTests.hpp ErrorHandlingForTests.hpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BasicIdentifierAdapter
  SOURCES BasicIdentifierAdapter.cpp AdaptersForTests.hpp ErrorHandlingForTests.hpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )
