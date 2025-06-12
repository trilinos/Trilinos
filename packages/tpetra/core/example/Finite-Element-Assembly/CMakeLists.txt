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

set(COUNTER 1)
set(MAX_COUNT 11)
WHILE(COUNTER LESS MAX_COUNT)
  math(EXPR NP "${COUNTER} * ${COUNTER}")
  TRIBITS_ADD_TEST(
      FEMAssembly
      NAME Performance_StrongScaling_FEMAssembly_InsertGlobalIndicesFESPKokkos
      ARGS "--with-insert-global-indices-fe --num-elements-x=8192 --num-elements-y=8192 --kokkos"
      COMM mpi
      NUM_MPI_PROCS ${NP}
      STANDARD_PASS_OUTPUT
      RUN_SERIAL
      CATEGORIES PERFORMANCE
  )
  math(EXPR COUNTER "${COUNTER} + 1")
ENDWHILE()

