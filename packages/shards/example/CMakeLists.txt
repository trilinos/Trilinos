
INCLUDE(PackageAddExecutableAndTest)

IF(Shards_ENABLE_Teuchos)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    example_01
    SOURCES example_01.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

ENDIF()

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  example_02
  SOURCES example_02.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  example_03
  SOURCES example_03.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
