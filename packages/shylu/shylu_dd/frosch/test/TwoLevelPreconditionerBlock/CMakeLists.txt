TRIBITS_ADD_EXECUTABLE(
    twolevelblock
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(TwoLevelBlockCopyFiles
    DEST_FILES GDSW.xml RGDSW.xml IPOUHarmonic.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS twolevelblock
)

###########################
# 2D - 1 DPN - Ordering 0 #
###########################
TRIBITS_ADD_TEST(
    twolevelblock
    NAME test_twolevelblock_2d_1dpn_o0_gdsw
    ARGS "--Dim=2 --DPN=1 --Ordering=0 --List=GDSW.xml"
    COMM mpi
    NUM_MPI_PROCS 4
)