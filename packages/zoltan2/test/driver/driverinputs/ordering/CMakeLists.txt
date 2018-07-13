IF (${PACKAGE_NAME}_ENABLE_Scotch)
TRIBITS_ADD_TEST(
    test_driver
    NAME orderingTestDriverExample
    DIRECTORY ../..
    NUM_MPI_PROCS 1
    COMM serial mpi
    ARGS
    "orderingTest.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_ordering
    SOURCE_FILES
        orderingTest.mtx
        orderingTest.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)