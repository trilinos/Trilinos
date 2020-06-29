# Test problem that exposes a divide-by-zeroin MJ when num_global_parts==16
TRIBITS_ADD_TEST(
    test_driver
    NAME sacer4320
    DIRECTORY ../..
    NUM_MPI_PROCS 8
    COMM mpi
    ARGS
    "sacer.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_sacer4320
    SOURCE_FILES
        sacer.xml
        sacer.graph
        sacer.assign
        sacer.coords
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
