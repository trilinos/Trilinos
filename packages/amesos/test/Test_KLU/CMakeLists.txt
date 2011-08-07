
INCLUDE(PackageAddExecutableAndTest)


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Galeri)
IF (${PACKAGE_NAME}_ENABLE_Galeri)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Test_KLU
    SOURCES cxx_main.cpp
    ARGS 
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

  # Added to excercise a memory leak in 
  # Amesos_KLU associated with reindexing
  # test is retained for primarily historical
  # reasons
  #PACKAGE_ADD_EXECUTABLE_AND_TEST(
  #  Test_KLU_memleak
  #  SOURCES cxx_memleak_main.cpp
  #  ARGS 
  #  NUM_MPI_PROCS 4
  #  COMM serial mpi
  #  )

ENDIF()
