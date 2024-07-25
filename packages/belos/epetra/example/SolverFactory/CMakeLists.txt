
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Galeri)

IF (${PACKAGE_NAME}_ENABLE_Galeri)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SolverFactory_Epetra_Galeri_Ex
    SOURCES SolverFactoryEpetraGaleriEx.cpp
    COMM serial mpi
    ARGS --verbose
    STANDARD_PASS_OUTPUT
    )

ENDIF (${PACKAGE_NAME}_ENABLE_Galeri)

