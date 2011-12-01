
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../Utilities)

ASSERT_DEFINED(Phalanx_ENABLE_Intrepid)

IF (Phalanx_ENABLE_Intrepid)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    intrepid_interoperability
    SOURCES IntrepidInteroperability.cpp
    DEPLIBS phalanx_test_utilities
    COMM serial mpi
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

ENDIF(Phalanx_ENABLE_Intrepid)