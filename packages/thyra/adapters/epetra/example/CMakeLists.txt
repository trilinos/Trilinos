INCLUDE(PackageAddExecutableAndTest)


SET(COMMON_SOURCES createTridiagEpetraLinearOp.cpp)

INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PARENT_PACKAGE_SOURCE_DIR}/core/example/operator_vector)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  sillyCgSolve_epetra
  SOURCES sillyCgSolve_epetra.cpp ${COMMON_SOURCES}
  COMM serial mpi
  ARGS --verbose
  )



PACKAGE_ADD_EXECUTABLE_AND_TEST(
  sillyPowerMethod_epetra
  SOURCES sillyPowerMethod_epetra.cpp ${COMMON_SOURCES}
  COMM serial mpi
 )
