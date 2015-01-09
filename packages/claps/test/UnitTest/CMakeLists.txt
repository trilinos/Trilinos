

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

IF (TPL_ENABLE_MPI)
  IF (HAVE_SCLAPS)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
        Instantiation_UnitTest
        SOURCES Claps_Instantiation_UnitTest.cpp Claps_UnitTest.cpp 
        DEPLIBS claps
        NUM_MPI_PROCS 1
        STANDARD_PASS_OUTPUT
        )
  ENDIF()
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TrivialSerial_UnitTest
    SOURCES Claps_TrivialSerial_UnitTest.cpp Claps_UnitTest.cpp 
    DEPLIBS claps
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )


