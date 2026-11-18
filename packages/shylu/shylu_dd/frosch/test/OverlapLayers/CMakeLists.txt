TRIBITS_ADD_EXECUTABLE(
overlap_layers
SOURCES main.cpp
)

################################################
# ExtendOverlapByLayers_Threaded vs references #
################################################
TRIBITS_ADD_TEST(
overlap_layers
NAME test_overlap_layers_chain_NP1
ARGS "--case=chain"
COMM mpi
NUM_MPI_PROCS 1
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_layers
NAME test_overlap_layers_chain_NP2
ARGS "--case=chain"
COMM mpi
NUM_MPI_PROCS 2
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_layers
NAME test_overlap_layers_chain_NP4
ARGS "--case=chain"
COMM mpi
NUM_MPI_PROCS 4
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_layers
NAME test_overlap_layers_grid_NP4
ARGS "--case=grid"
COMM mpi
NUM_MPI_PROCS 4
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_layers
NAME test_overlap_layers_empty_NP3
ARGS "--case=empty"
COMM mpi
NUM_MPI_PROCS 3
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_layers
NAME test_overlap_layers_disconnected_NP3
ARGS "--case=disconnected"
COMM mpi
NUM_MPI_PROCS 3
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_layers
NAME test_overlap_layers_nodiag_NP2
ARGS "--case=nodiag"
COMM mpi
NUM_MPI_PROCS 2
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_layers
NAME test_overlap_layers_gapped_NP4
ARGS "--case=gapped"
COMM mpi
NUM_MPI_PROCS 4
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_layers
NAME test_overlap_layers_rehash_NP4
ARGS "--case=rehash"
COMM mpi
NUM_MPI_PROCS 4
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_layers
NAME test_overlap_layers_galeri_DIM2_DPN1_NP4
ARGS "--case=galeri --DIM=2 --DPN=1 --M=5"
COMM mpi
NUM_MPI_PROCS 4
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_layers
NAME test_overlap_layers_galeri_DIM2_DPN2_NP4
ARGS "--case=galeri --DIM=2 --DPN=2 --M=5"
COMM mpi
NUM_MPI_PROCS 4
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_layers
NAME test_overlap_layers_galeri_DIM3_DPN1_NP8
ARGS "--case=galeri --DIM=3 --DPN=1 --M=5"
COMM mpi
NUM_MPI_PROCS 8
ENVIRONMENT OMP_NUM_THREADS=4
)

TRIBITS_ADD_TEST(
overlap_layers
NAME test_overlap_layers_galeri_DIM3_DPN3_NP8
ARGS "--case=galeri --DIM=3 --DPN=3 --M=5"
COMM mpi
NUM_MPI_PROCS 8
ENVIRONMENT OMP_NUM_THREADS=4
)
