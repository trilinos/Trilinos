TRIBITS_ADD_EXECUTABLE(
    twolevelinputfiles
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(TwoLevelInputFilesCopyFiles
    DEST_FILES Parameters.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS twolevelinputfiles
)

#TRIBITS_ADD_TEST(
#    twolevelinputfiles
#    NAME test_twolevelinputfiles_2d_gdsw
#    ARGS "--DIM=2 --CL=GDSW"
#    COMM mpi
#    NUM_MPI_PROCS 4
#)
