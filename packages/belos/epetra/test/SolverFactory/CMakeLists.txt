
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Galeri)
IF (${PACKAGE_NAME}_ENABLE_Galeri)
    
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    solver_factory_test
    SOURCES test_solver_factory.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

ENDIF (${PACKAGE_NAME}_ENABLE_Galeri)
