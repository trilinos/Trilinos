
TRIBITS_ADD_EXECUTABLE_AND_TEST(
    FECrsGraph_UnitTests
    SOURCES
        FECrsGraph_UnitTests.cpp
        ${TEUCHOS_STD_UNIT_TEST_MAIN}
    ARGS ""
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
