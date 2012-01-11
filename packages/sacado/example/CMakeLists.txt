INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  dfad_example
  SOURCES dfad_example.cpp
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ad_example
  SOURCES ad_example.cpp
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

IF(Sacado_ENABLE_NESTEDFAD)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    dfad_dfad_example
    SOURCES dfad_dfad_example.cpp
    ARGS -v
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  sfad_example
  SOURCES sfad_example.cpp
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  dmfad_example
  SOURCES dmfad_example.cpp
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  trad_example
  SOURCES trad_example.cpp
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tradvec_example
  SOURCES tradvec_example.cpp
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  trad_dfad_example
  SOURCES trad_dfad_example.cpp
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  taylor_example
  SOURCES taylor_example.cpp
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  dfad_taylor_example
  SOURCES dfad_taylor_example.cpp
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  dfad_sfc_example
  SOURCES dfad_sfc_example.cpp
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  trad_sfc_example
  SOURCES trad_sfc_example.cpp
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  sparse_example
  SOURCES sparse_example.cpp
  ARGS -v
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

IF(Sacado_ENABLE_Teuchos)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    blas_example
    SOURCES blas_example.cpp
    ARGS -v
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    vector_blas_example
    SOURCES vector_blas_example.cpp
    ARGS -v
    COMM serial mpi
    NUM_MPI_PROCS 1
    )
ENDIF()
