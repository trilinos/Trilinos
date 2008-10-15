INCLUDE(Trilinos_Add_Executable_And_Test)


SET(COMMON_SOURCES createTridiagEpetraLinearOp.cpp)



TRILINOS_ADD_EXECUTABLE_AND_TEST(
  sillyCgSolve_epetra
  SOURCES sillyCgSolve_epetra.cpp ${COMMON_SOURCES}
  COMM serial mpi
  ARGS --verbose
  )



TRILINOS_ADD_EXECUTABLE_AND_TEST(
  sillyPowerMethod_epetra
  SOURCES sillyPowerMethod_epetra.cpp ${COMMON_SOURCES}
  COMM serial mpi
  ARGS --verbose 
 )
