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
  ARGS "--with-insert-global-indices-fe --num-elements-x=200 --num-elements-y=200 --with-kokkos --with-compare-graph"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)

# The Kokkos-assembled graph must be identical to the inserted graph on a
# single process too (no owned+shared -> owned communication).
TRIBITS_ADD_TEST(
  FEMAssembly
  NAME FEMAssembly_InsertGlobalIndicesFESPKokkosOneRank
  ARGS "--with-insert-global-indices-fe --num-elements-x=50 --num-elements-y=50 --with-kokkos --with-compare-graph"
  COMM serial mpi
  NUM_MPI_PROCS 1
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
  ARGS "--with-total-element-loop --num-elements-x=200 --num-elements-y=200 --with-kokkos"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)

###################################
# Performance Testing
###################################

SET(COUNTER 1)
SET(MAX_COUNT 11)
WHILE(COUNTER LESS MAX_COUNT)
  MATH(EXPR NP "${COUNTER} * ${COUNTER}")
  TRIBITS_ADD_TEST(
    FEMAssembly
    NAME Performance_StrongScaling_FEMAssembly_InsertGlobalIndicesFESPKokkos
    ARGS "--with-insert-global-indices-fe --num-elements-x=8192 --num-elements-y=8192 --with-kokkos"
    COMM mpi
    NUM_MPI_PROCS ${NP}
    STANDARD_PASS_OUTPUT
    RUN_SERIAL
    CATEGORIES PERFORMANCE
  )
  MATH(EXPR COUNTER "${COUNTER} + 1")
ENDWHILE()
