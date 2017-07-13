TRIBITS_ADD_TEST(
    test_driver
    NAME rectilinearRCB_left
    DIRECTORY ../..
    NUM_MPI_PROCS 2
    COMM serial mpi
    ARGS
    "rectilinearRCB_left.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    test_driver
    NAME rectilinearRCB_right
    DIRECTORY ../..
    NUM_MPI_PROCS 2
    COMM serial mpi
    ARGS
    "rectilinearRCB_right.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_rectilinearRCB
    SOURCE_FILES
        rectilinearRCB_left.xml
        rectilinearRCB_left.graph
        rectilinearRCB_left.coords
        rectilinearRCB_right.xml
        rectilinearRCB_right.graph
        rectilinearRCB_right.coords
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
