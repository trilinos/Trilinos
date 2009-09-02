
INCLUDE(PackageAddExecutableAndTest)


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_SuperLUDist)
IF (${PACKAGE_NAME}_ENABLE_SuperLUDist)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Test_SuperLU_Dist
    SOURCES cxx_main.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

ENDIF()
