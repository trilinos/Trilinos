IF(PANZER_HAVE_EPETRA_STACK)

  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  SET(UNIT_TEST_DRIVER
    ${PANZER_UNIT_TEST_MAIN})

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    periodic_bcs
    SOURCES periodic_bcs.cpp ${UNIT_TEST_DRIVER}
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    periodic_mesh
    SOURCES periodic_mesh.cpp ${UNIT_TEST_DRIVER}
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  IF(PANZER_HAVE_STKSEARCH)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      periodic_search
      SOURCES periodic_search.cpp ${UNIT_TEST_DRIVER}
      COMM serial mpi
      NUM_MPI_PROCS 2
      )
  ENDIF()

  # TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #   periodic_mesh_nosubcells
  #   SOURCES periodic_mesh_nosubcells.cpp ${UNIT_TEST_DRIVER}
  #   COMM serial mpi
  #   NUM_MPI_PROCS 1
  #   )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    periodic_32bit_int_limit
    SOURCES periodic_32bit_int_limit.cpp ${UNIT_TEST_DRIVER}
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(periodic_32bit_int_limit_files
    SOURCE_FILES periodic_32bit_int_limit.jou periodic_32bit_int_limit.exo README.txt
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/"
    EXEDEPS periodic_32bit_int_limit
    )

ENDIF(PANZER_HAVE_EPETRA_STACK)
