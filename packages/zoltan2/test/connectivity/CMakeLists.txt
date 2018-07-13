
TRIBITS_ADD_EXECUTABLE_AND_TEST(
    iceSheet
    SOURCES iceSheetMain.cpp 
    NUM_MPI_PROCS 4
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

## MultiJagged Tests
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    mjTest
#    SOURCES MultiJaggedTest.cpp
#    NUM_MPI_PROCS 4
#    COMM serial mpi
#    ARGS
#    "I=1.03 P=2,2 O=0 C=4 F=simple TB=1"
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#)


#TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_more_files_for_partition_tests
#  SOURCE_FILES
#    Poisson.xml
#    torus.dmg
#    4/0.smb
#    4/1.smb	
#    4/2.smb
#    4/3.smb
#    pumiTri14/plate.dmg
#    pumiTri14/plate0.smb
#    pumiTri14/2/0.smb
#    pumiTri14/2/1.smb
#  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
#  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
#)
