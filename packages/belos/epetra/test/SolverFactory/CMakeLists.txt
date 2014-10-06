
ASSERT_DEFINED(Belos_ENABLE_Galeri)
IF (Belos_ENABLE_Galeri)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    solver_factory_test
    SOURCES test_solver_factory.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

ENDIF (Belos_ENABLE_Galeri)
