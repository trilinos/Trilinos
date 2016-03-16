TRIBITS_ADD_TEST(
    test_driver
    NAME glaze2RCB
    DIRECTORY ../..
    NUM_MPI_PROCS 2
    COMM serial mpi
    ARGS
    "glaze2.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_glaze2
    SOURCE_FILES
        glaze2.xml
        glaze2.graph
        glaze2.coords
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
