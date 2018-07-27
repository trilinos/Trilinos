TRIBITS_ADD_EXECUTABLE(
   iceSheet
   SOURCES ${HEADERS} iceSheetMain.cpp
   COMM serial mpi
)


TRIBITS_ADD_TEST(
    iceSheet
    NAME iceSheet_Antarctica_1K
    ARGS "meshfilename.txt" 
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    iceSheet
    NAME iceSheet_Checkerboard
    ARGS "checkboard.txt" 
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)


#TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_more_files_for_partition_tests
#  SOURCE_FILES
#    checkboard.txt
#    meshfilename.txt 
#  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/meshfiles
#  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
#)
