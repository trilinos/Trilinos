BELOS_ADD_EXECUTABLE_AND_TEUCHOS_AND_KOKKOS_TESTS(
  pcpg_tpetraex
  SOURCES test_pcpg_tpetraex.cpp
  COMM mpi
  NUM_MPI_PROCS 4
  )
