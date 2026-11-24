TRIBITS_ADD_EXECUTABLE(
overlap_strategies
SOURCES main.cpp
)

#################################################
# Adding Layers Strategy: Threaded vs CrsGraph  #
#################################################
TRIBITS_ADD_TEST(
overlap_strategies
NAME test_overlap_strategies_DIM2_DPN1_O0_TPETRA
ARGS "--M=5 --DIM=2 --O=0 --DPN=1"
COMM mpi
NUM_MPI_PROCS 4
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_strategies
NAME test_overlap_strategies_DIM3_DPN1_O0_TPETRA
ARGS "--M=5 --DIM=3 --O=0 --DPN=1"
COMM mpi
NUM_MPI_PROCS 8
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_strategies
NAME test_overlap_strategies_DIM2_DPN1_O1_TPETRA
ARGS "--M=5 --DIM=2 --O=1 --DPN=1"
COMM mpi
NUM_MPI_PROCS 4
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_strategies
NAME test_overlap_strategies_DIM2_DPN1_O2_TPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=1"
COMM mpi
NUM_MPI_PROCS 4
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_strategies
NAME test_overlap_strategies_DIM2_DPN2_O1_TPETRA
ARGS "--M=5 --DIM=2 --O=1 --DPN=2"
COMM mpi
NUM_MPI_PROCS 4
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_strategies
NAME test_overlap_strategies_DIM2_DPN2_O2_TPETRA
ARGS "--M=5 --DIM=2 --O=2 --DPN=2"
COMM mpi
NUM_MPI_PROCS 4
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_strategies
NAME test_overlap_strategies_DIM3_DPN1_O1_TPETRA
ARGS "--M=5 --DIM=3 --O=1 --DPN=1"
COMM mpi
NUM_MPI_PROCS 8
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_strategies
NAME test_overlap_strategies_DIM3_DPN1_O2_TPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=1"
COMM mpi
NUM_MPI_PROCS 8
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_strategies
NAME test_overlap_strategies_DIM3_DPN3_O1_TPETRA
ARGS "--M=5 --DIM=3 --O=1 --DPN=3"
COMM mpi
NUM_MPI_PROCS 8
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_strategies
NAME test_overlap_strategies_DIM3_DPN3_O2_TPETRA
ARGS "--M=5 --DIM=3 --O=2 --DPN=3"
COMM mpi
NUM_MPI_PROCS 8
ENVIRONMENT OMP_NUM_THREADS=4
)
