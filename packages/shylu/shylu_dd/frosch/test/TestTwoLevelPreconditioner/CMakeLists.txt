TRIBITS_ADD_EXECUTABLE(
    twolevel
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(TwoLevelCopyFiles
    DEST_FILES Parameters.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS twolevel
)

TRIBITS_ADD_TEST(
    twolevel
    NAME test_twolevel_2d_gdsw
    ARGS "--DIM=2 --CL=GDSW"
    COMM mpi
    NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
    twolevel
    NAME test_twolevel_2d_rgdsw
    ARGS "--DIM=2 --CL=RGDSW"
    COMM mpi
    NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
    twolevel
    NAME test_twolevel_2d_ipouharmonic
    ARGS "--DIM=2 --CL=IPOUHarmonic"
    COMM mpi
    NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
    twolevel
    NAME test_twolevel_3d_gdsw
    ARGS "--DIM=3 --CL=GDSW"
    COMM mpi
    NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
    twolevel
    NAME test_twolevel_3d_rgdsw
    ARGS "--DIM=3 --CL=RGDSW"
    COMM mpi
    NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
    twolevel
    NAME test_twolevel_3d_ipouharmonic
    ARGS "--DIM=3 --CL=IPOUHarmonic"
    COMM mpi
    NUM_MPI_PROCS 8
)
