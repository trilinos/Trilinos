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
    NAME iceSheet_test1
    ARGS " test1.quad.msh test1-borders.quad.msh test1-basal-friction.ascii "
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)


TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_files_for_connectivity_tests
  SOURCE_FILES
    test1.quad.msh
    test1-borders.quad.msh
    test1-basal-friction.ascii
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/meshData
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
