

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    pcpg_tpetraex
    SOURCES test_pcpg_tpetraex.cpp
    COMM mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    pcpg_densekokkos_tpetraex
    SOURCES test_kokkos_pcpg_tpetraex.cpp
    COMM mpi
    NUM_MPI_PROCS 4
    )

ENDIF()


