TRIBITS_ADD_EXECUTABLE(
    twolevelblockinputfiles
    SOURCES main.cpp
    
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(TwoLevelBlockInputFilesCopyFiles
    DEST_FILES Parameters.xml stokes.h5
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS twolevelblockinputfiles
)

TRIBITS_ADD_TEST(
    twolevelblockinputfiles
    NAME test_twolevelblockinputfiles_2d_gdsw
    ARGS
    COMM mpi
    NUM_MPI_PROCS 4
)
