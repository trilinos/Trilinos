# Test problem that exposes a memory leak in MJ when num_global_parts==1
TRIBITS_ADD_TEST(
    test_driver
    NAME coffeyMJ
    DIRECTORY ../..
    NUM_MPI_PROCS 1
    COMM serial mpi
    ARGS
    "coffey.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_coffey
    SOURCE_FILES
        coffey.xml
        coffey.graph
        coffey.coords
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
