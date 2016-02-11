TRIBITS_ADD_TEST(
    test_driver
    NAME rectilinearRCB
    DIRECTORY ../..
    NUM_MPI_PROCS 2
    COMM serial mpi
    ARGS
    "rectilinearRCB.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_rectilinearRCB
    SOURCE_FILES
        rectilinearRCB.xml
        rectilinearRCB.graph
        rectilinearRCB.coords
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
