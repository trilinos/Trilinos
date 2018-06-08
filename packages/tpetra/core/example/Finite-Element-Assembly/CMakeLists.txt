
# This example works whether or not MPI is enabled.
# It does not refer to MPI explicitly.
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  FEMAssembly
#  SOURCES fem_assembly_main.cpp
#  ARGS 
#  POSTFIX_AND_ARGS_0 InsertGlobalIndicesDP
#  "--with-insert-global-indices-dp --num-elements-x=1000 --num-elements-y=1000"
#  POSTFIX_AND_ARGS_1 LocalElementLoopDP
#  "--with-local-element-loop-dp --num-elements-x=1000 --num-elements-y=1000"
#  POSTFIX_AND_ARGS_2 TotalElementLoopDP
#  "--with-total-element-loop-dp --num-elements-x=1000 --num-elements-y=1000"
#  POSTFIX_AND_ARGS_3 TotalElementLoopSP
#  "--with-total-element-loop-sp --num-elements-x=1000 --num-elements-y=1000"
#  COMM  mpi
#  NUM_MPI_PROCS 4
#  STANDARD_PASS_OUTPUT
#)

TRIBITS_ADD_EXECUTABLE(
  FEMAssembly
  SOURCES fem_assembly_main.cpp
  COMM serial mpi
)

TRIBITS_ADD_TEST(
  FEMAssembly
  NAME FEMAssembly_InsertGlobalIndicesDP
  ARGS "--with-insert-global-indices --with-DynamicProfile --num-elements-x=200 --num-elements-y=200"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  FEMAssembly
  NAME FEMAssembly_LocalElementLoopDP
  ARGS "--with-local-element-loop --with-DynamicProfile --num-elements-x=200 --num-elements-y=200"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  FEMAssembly
  NAME FEMAssembly_TotalElementLoopDP
  ARGS "--with-total-element-loop --with-DynamicProfile --num-elements-x=200 --num-elements-y=200"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  FEMAssembly
  NAME FEMAssembly_TotalElementLoopSP
  ARGS "--with-total-element-loop --with-StaticProfile --num-elements-x=200 --num-elements-y=200"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)

