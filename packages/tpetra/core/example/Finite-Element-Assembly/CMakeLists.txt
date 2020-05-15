TRIBITS_ADD_EXECUTABLE(
  FEMAssembly
  SOURCES fem_assembly_main.cpp
  COMM serial mpi
  CATEGORIES BASIC PERFORMANCE
)

TRIBITS_ADD_TEST(
    FEMAssembly
    NAME FEMAssembly_InsertGlobalIndicesFESP
    ARGS "--with-insert-global-indices-fe --num-elements-x=200 --num-elements-y=200"
    COMM serial mpi
    NUM_MPI_PROCS 4
    STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
    FEMAssembly 
    NAME FEMAssembly_InsertGlobalIndicesFESPKokkos
    ARGS "--with-insert-global-indices-fe --num-elements-x=200 --num-elements-y=200 --kokkos"
    COMM serial mpi
    NUM_MPI_PROCS 4
    STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  FEMAssembly
  NAME FEMAssembly_TotalElementLoopSP
  ARGS "--with-total-element-loop --num-elements-x=200 --num-elements-y=200"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  FEMAssembly
  NAME FEMAssembly_TotalElementLoopSPKokkos
  ARGS "--with-total-element-loop --num-elements-x=200 --num-elements-y=200 --kokkos"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)

###################################
# Performance Testing
###################################

TRIBITS_ADD_TEST(
  FEMAssembly
  NAME Performance_StrongScaling_FEMAssembly_InsertGlobalIndicesFESPKokkos
  ARGS "--with-insert-global-indices-fe --num-elements-x=8192 --num-elements-y=8192 --kokkos"
  COMM mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  RUN_SERIAL
  CATEGORIES PERFORMANCE
)

TRIBITS_ADD_TEST(
  FEMAssembly
  NAME Performance_StrongScaling_FEMAssembly_InsertGlobalIndicesFESPKokkos
  ARGS "--with-insert-global-indices-fe --num-elements-x=8192 --num-elements-y=8192 --kokkos"
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  RUN_SERIAL
  CATEGORIES PERFORMANCE
)

TRIBITS_ADD_TEST(
  FEMAssembly
  NAME Performance_StrongScaling_FEMAssembly_InsertGlobalIndicesFESPKokkos
  ARGS "--with-insert-global-indices-fe --num-elements-x=8192 --num-elements-y=8192 --kokkos"
  COMM mpi
  NUM_MPI_PROCS 9
  STANDARD_PASS_OUTPUT
  RUN_SERIAL
  CATEGORIES PERFORMANCE
)

TRIBITS_ADD_TEST(
  FEMAssembly
  NAME Performance_StrongScaling_FEMAssembly_InsertGlobalIndicesFESPKokkos
  ARGS "--with-insert-global-indices-fe --num-elements-x=8192 --num-elements-y=8192 --kokkos"
  COMM mpi
  NUM_MPI_PROCS 16
  STANDARD_PASS_OUTPUT
  RUN_SERIAL
  CATEGORIES PERFORMANCE
)

TRIBITS_ADD_TEST(
  FEMAssembly
  NAME Performance_StrongScaling_FEMAssembly_InsertGlobalIndicesFESPKokkos
  ARGS "--with-insert-global-indices-fe --num-elements-x=8192 --num-elements-y=8192 --kokkos"
  COMM mpi
  NUM_MPI_PROCS 25
  STANDARD_PASS_OUTPUT
  RUN_SERIAL
  CATEGORIES PERFORMANCE
)

