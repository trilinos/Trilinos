
INCLUDE(TribitsAddExecutableAndTest)


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_AztecOO)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Galeri)

IF (NOT Trilinos_NO_32BIT_GLOBAL_INDICES)

IF (${PACKAGE_NAME}_ENABLE_Galeri)
  IF (${PACKAGE_NAME}_ENABLE_AztecOO)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      ex_Factory
      SOURCES Ifpack_ex_Factory.cpp
      ARGS "-v"
      NUM_MPI_PROCS 1-4
      COMM serial mpi
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      ex_ICT
      SOURCES Ifpack_ex_ICT.cpp
      ARGS "-v"
      NUM_MPI_PROCS 1-4
      COMM serial mpi
      )

    ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos)
    IF (${PACKAGE_NAME}_ENABLE_Amesos)

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        ex_Amesos
        SOURCES Ifpack_ex_Amesos.cpp
        ARGS "-v"
        NUM_MPI_PROCS 1-4
        COMM serial mpi
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        ex_BlockRelaxation
        SOURCES Ifpack_ex_BlockRelaxation.cpp
        ARGS "-v"
        NUM_MPI_PROCS 1-4
        COMM serial mpi
        )

    ENDIF()

  ENDIF()

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ex_reordering
    SOURCES Ifpack_ex_Reordering.cpp
    ARGS "-v"
    NUM_MPI_PROCS 1-4
    COMM serial mpi
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ex_filtering
    SOURCES Ifpack_ex_Filtering.cpp
    ARGS "-v"
    NUM_MPI_PROCS 1-4
    COMM serial mpi
    )

ENDIF()

ENDIF()

IF (NOT Trilinos_NO_64BIT_GLOBAL_INDICES)

IF (${PACKAGE_NAME}_ENABLE_Galeri)
  IF (${PACKAGE_NAME}_ENABLE_AztecOO)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      ex_Factory_LL
      SOURCES Ifpack_ex_Factory_LL.cpp
      ARGS "-v"
      NUM_MPI_PROCS 1-4
      COMM serial mpi
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      ex_ICT_LL
      SOURCES Ifpack_ex_ICT_LL.cpp
      ARGS "-v"
      NUM_MPI_PROCS 1-4
      COMM serial mpi
      )

    ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos)
    IF (${PACKAGE_NAME}_ENABLE_Amesos)

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        ex_Amesos_LL
        SOURCES Ifpack_ex_Amesos_LL.cpp
        ARGS "-v"
        NUM_MPI_PROCS 1-4
        COMM serial mpi
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        ex_BlockRelaxation_LL
        SOURCES Ifpack_ex_BlockRelaxation_LL.cpp
        ARGS "-v"
        NUM_MPI_PROCS 1-4
        COMM serial mpi
        )

    ENDIF()

  ENDIF()

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ex_reordering_LL
    SOURCES Ifpack_ex_Reordering_LL.cpp
    ARGS "-v"
    NUM_MPI_PROCS 1-4
    COMM serial mpi
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ex_filtering_LL
    SOURCES Ifpack_ex_Filtering_LL.cpp
    ARGS "-v"
    NUM_MPI_PROCS 1-4
    COMM serial mpi
    )

ENDIF()

ENDIF()
