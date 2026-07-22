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
## Tpetra
TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Graph_DIM2_DPN1_TPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Graph.xml"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Graph_DIM2_DPN2_TPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=2 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Graph.xml"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Graph_DIM3_DPN1_TPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Graph.xml"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Graph_DIM3_DPN3_TPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=3 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Graph.xml"
COMM mpi
NUM_MPI_PROCS 8
)

################################
# Overlap Determination Matrix #
################################
## Tpetra
TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Matrix_DIM2_DPN1_TPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Matrix.xml"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Matrix_DIM2_DPN2_TPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=2 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Matrix.xml"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Matrix_DIM3_DPN1_TPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Matrix.xml"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Matrix_DIM3_DPN3_TPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=3 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Matrix.xml"
COMM mpi
NUM_MPI_PROCS 8
)

#############################
# Overlap Determination Old #
#############################
## Tpetra
TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Old_DIM2_DPN1_TPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Old.xml"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Old_DIM2_DPN2_TPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=2 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Old.xml"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Old_DIM3_DPN1_TPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=1 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Old.xml"
COMM mpi
NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
overlap
NAME test_overlap_Old_DIM3_DPN3_TPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=3 --PLIST=ParameterLists/ParameterList_OneLevelPreconditioner_Old.xml"
COMM mpi
NUM_MPI_PROCS 8
)
