INCLUDE(TribitsAddExecutableAndTest)


SET(COMMON_SOURCES createTridiagEpetraLinearOp.cpp)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PARENT_PACKAGE_SOURCE_DIR}/core/example/operator_vector)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  sillyCgSolve_epetra
  SOURCES sillyCgSolve_epetra.cpp ${COMMON_SOURCES}
  COMM serial mpi
  ARGS --verbose
  )



TRIBITS_ADD_EXECUTABLE_AND_TEST(
  sillyPowerMethod_epetra
  SOURCES sillyPowerMethod_epetra.cpp ${COMMON_SOURCES}
  COMM serial mpi
 )
