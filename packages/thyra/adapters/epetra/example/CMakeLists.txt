INCLUDE(PackageAddExecutableAndTest)


SET(COMMON_SOURCES createTridiagEpetraLinearOp.cpp)



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
  ARGS --verbose 
 )
