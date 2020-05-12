# Test problem that exposes a memory leak in MJ when num_global_parts==1
TRIBITS_ADD_TEST(
    test_driver
    NAME kuberry
    DIRECTORY ../..
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "kuberry.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_kuberry
    SOURCE_FILES
        kuberry.xml
        kuberry.graph
        kuberry.coords
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
