TRIBITS_ADD_EXECUTABLE(
    twolevel
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(TwoLevelCopyFiles
    DEST_FILES IPOUHarmonic.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS twolevel
)

###########################
# 2D - 1 DPN - Ordering 0 #
###########################
#TRIBITS_ADD_TEST(
#    twolevel
#    NAME test_twolevel_2d_1dpn_o0_gdsw
#    ARGS "--Dim=2 --DPN=1 --Ordering=0 --List=GDSW.xml"
#    COMM mpi
#    NUM_MPI_PROCS 4
#)

#TRIBITS_ADD_TEST(
#    twolevel
#    NAME test_twolevel_2d_1dpn_o0_rgdsw
#    ARGS "--Dim=2 --DPN=1 --Ordering=0 --List=RGDSW.xml"
#    COMM mpi
#    NUM_MPI_PROCS 4
#)

TRIBITS_ADD_TEST(
    twolevel
    NAME test_twolevel_2d_1dpn_o0_ipouharmonic
    ARGS "--Dim=2 --DPN=1 --Ordering=0 --List=IPOUHarmonic.xml"
    COMM mpi
    NUM_MPI_PROCS 4
)

###########################
# 2D - 1 DPN - Ordering 1 #
###########################
#TRIBITS_ADD_TEST(
#    twolevel
#    NAME test_twolevel_2d_1dpn_o1_gdsw
#    ARGS "--Dim=2 --DPN=1 --Ordering=0 --List=GDSW.xml"
#    COMM mpi
#    NUM_MPI_PROCS 4
#)

#TRIBITS_ADD_TEST(
#    twolevel
#    NAME test_twolevel_2d_1dpn_o1_rgdsw
#    ARGS "--Dim=2 --DPN=1 --Ordering=0 --List=RGDSW.xml"
#    COMM mpi
#    NUM_MPI_PROCS 4
#)

TRIBITS_ADD_TEST(
    twolevel
    NAME test_twolevel_2d_1dpn_o1_ipouharmonic
    ARGS "--Dim=2 --DPN=1 --Ordering=0 --List=IPOUHarmonic.xml"
    COMM mpi
    NUM_MPI_PROCS 4
)

###########################
# 2D - 2 DPN - Ordering 0 #
###########################
#TRIBITS_ADD_TEST(
#    twolevel
#    NAME test_twolevel_2d_2dpn_o0_gdsw
#    ARGS "--Dim=2 --DPN=2 --Ordering=0 --List=GDSW.xml"
#    COMM mpi
#    NUM_MPI_PROCS 4
#)

#TRIBITS_ADD_TEST(
#    twolevel
#    NAME test_twolevel_2d_2dpn_o0_rgdsw
#    ARGS "--Dim=2 --DPN=2 --Ordering=0 --List=RGDSW.xml"
#    COMM mpi
#    NUM_MPI_PROCS 4
#)

TRIBITS_ADD_TEST(
    twolevel
    NAME test_twolevel_2d_2dpn_o0_ipouharmonic
    ARGS "--Dim=2 --DPN=2 --Ordering=0 --List=IPOUHarmonic.xml"
    COMM mpi
    NUM_MPI_PROCS 4
)

###########################
# 2D - 2 DPN - Ordering 1 #
###########################
#TRIBITS_ADD_TEST(
#    twolevel
#    NAME test_twolevel_2d_2dpn_o1_gdsw
#    ARGS "--Dim=2 --DPN=2 --Ordering=0 --List=GDSW.xml"
#    COMM mpi
#    NUM_MPI_PROCS 4
#)

#TRIBITS_ADD_TEST(
#    twolevel
#    NAME test_twolevel_2d_2dpn_o1_rgdsw
#    ARGS "--Dim=2 --DPN=2 --Ordering=0 --List=RGDSW.xml"
#    COMM mpi
#    NUM_MPI_PROCS 4
#)

TRIBITS_ADD_TEST(
    twolevel
    NAME test_twolevel_2d_2dpn_o1_ipouharmonic
    ARGS "--Dim=2 --DPN=2 --Ordering=0 --List=IPOUHarmonic.xml"
    COMM mpi
    NUM_MPI_PROCS 4
)
