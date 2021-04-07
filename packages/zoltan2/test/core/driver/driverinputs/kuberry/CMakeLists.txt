# Test problem for #7351 that exposes an out-of-bounds array access in MJ when rectilinear = 1
TRIBITS_ADD_TEST(
    test_driver
    NAME kuberry
    DIRECTORY ../..
    COMM serial mpi
    ARGS
    "kuberry.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_kuberry
    SOURCE_FILES
        kuberry.assign
        kuberry.xml
        kuberry.graph
        kuberry.coords
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
