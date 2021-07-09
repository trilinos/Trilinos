# Test problem that exposes a migration error in MJ when p=64
TRIBITS_ADD_TEST(
    test_driver
    NAME miniem4669
    DIRECTORY ../..
    NUM_MPI_PROCS 64 
    COMM mpi
    ARGS
    "miniem.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_miniem4669
    SOURCE_FILES
        miniem.xml
        miniem.graph
        miniem.coords
        miniem.assign
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
