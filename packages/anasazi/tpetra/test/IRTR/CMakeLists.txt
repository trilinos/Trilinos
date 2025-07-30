TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_IRTR_test
  SOURCES cxx_main.cpp
  ARGS
    "--sort=LR --verbose --skinny"
    "--sort=LR --verbose --hefty"
    "--sort=SR --verbose --skinny"
    "--sort=SR --verbose --hefty"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_IRTR_auxtest
  SOURCES cxx_main_auxtest.cpp
  ARGS
    "--sort=LR --verbose --skinny"
    "--sort=LR --verbose --hefty"
    "--sort=SR --verbose --skinny"
    "--sort=SR --verbose --hefty"
    "--verbose"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_IRTR_Lap_test
  SOURCES cxx_main_lap.cpp 
  ARGS 
  COMM serial mpi
  )

IF(Tpetra_INST_COMPLEX_DOUBLE)
ENDIF()
