
INCLUDE(PackageAddExecutableAndTest)


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_AztecOO)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Galeri)
IF (${PACKAGE_NAME}_ENABLE_Galeri)
  IF (${PACKAGE_NAME}_ENABLE_AztecOO)

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      ex_Factory
      SOURCES Ifpack_ex_Factory.cpp
      ARGS "-v"
      NUM_MPI_PROCS 1-4
      COMM serial mpi
      )

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      ex_ICT
      SOURCES Ifpack_ex_ICT.cpp
      ARGS "-v"
      NUM_MPI_PROCS 1-4
      COMM serial mpi
      )

    ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos)
    IF (${PACKAGE_NAME}_ENABLE_Amesos)

      PACKAGE_ADD_EXECUTABLE_AND_TEST(
        ex_Amesos
        SOURCES Ifpack_ex_Amesos.cpp
        ARGS "-v"
        NUM_MPI_PROCS 1-4
        COMM serial mpi
        )

      PACKAGE_ADD_EXECUTABLE_AND_TEST(
        ex_BlockRelaxation
        SOURCES Ifpack_ex_BlockRelaxation.cpp
        ARGS "-v"
        NUM_MPI_PROCS 1-4
        COMM serial mpi
        )

    ENDIF()

  ENDIF()

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    ex_reordering
    SOURCES Ifpack_ex_Reordering.cpp
    ARGS "-v"
    NUM_MPI_PROCS 1-4
    COMM serial mpi
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    ex_filtering
    SOURCES Ifpack_ex_Filtering.cpp
    ARGS "-v"
    NUM_MPI_PROCS 1-4
    COMM serial mpi
    )

ENDIF()
