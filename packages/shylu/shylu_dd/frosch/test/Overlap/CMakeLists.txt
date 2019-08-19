TRIBITS_ADD_EXECUTABLE(
overlap
SOURCES main.cpp
)

ADD_SUBDIRECTORIES(
ParameterLists
)

###############################
# Overlap Determination Graph #
###############################
## Epetra
IF(ShyLU_DDFROSch_ENABLE_Epetra AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Graph_DIM2_DPN1_EPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Graph.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Graph_DIM2_DPN2_EPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=2 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Graph.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Graph_DIM3_DPN1_EPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Graph.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Graph_DIM3_DPN3_EPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=3 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Graph.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Graph_DIM2_DPN1_TPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Graph.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Graph_DIM2_DPN2_TPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=2 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Graph.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Graph_DIM3_DPN1_TPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Graph.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Graph_DIM3_DPN3_TPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=3 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Graph.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

################################
# Overlap Determination Matrix #
################################
## Epetra
IF(ShyLU_DDFROSch_ENABLE_Epetra AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Matrix_DIM2_DPN1_EPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Matrix.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Matrix_DIM2_DPN2_EPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=2 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Matrix.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Matrix_DIM3_DPN1_EPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Matrix.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Matrix_DIM3_DPN3_EPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=3 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Matrix.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Matrix_DIM2_DPN1_TPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Matrix.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Matrix_DIM2_DPN2_TPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=2 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Matrix.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Matrix_DIM3_DPN1_TPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Matrix.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Matrix_DIM3_DPN3_TPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=3 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Matrix.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

#############################
# Overlap Determination Old #
#############################
## Epetra
IF(ShyLU_DDFROSch_ENABLE_Epetra AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Old_DIM2_DPN1_EPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Old.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Old_DIM2_DPN2_EPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=2 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Old.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Old_DIM3_DPN1_EPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Old.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Old_DIM3_DPN3_EPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=3 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Old.xml --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

## Tpetra
TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Old_DIM2_DPN1_TPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Old.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Old_DIM2_DPN2_TPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=2 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Old.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Old_DIM3_DPN1_TPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Old.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Old_DIM3_DPN3_TPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=3 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Old.xml --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
