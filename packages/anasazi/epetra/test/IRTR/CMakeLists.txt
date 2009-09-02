INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Epetra_IRTR_test
  SOURCES cxx_main.cpp
  ARGS 
    "--sort=LR --verbose --skinny"
    "--sort=LR --verbose --hefty"
    "--sort=SR --verbose --skinny"
    "--sort=SR --verbose --hefty"
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Epetra_IRTR_auxtest
  SOURCES cxx_main_auxtest.cpp
  ARGS 
    "--sort=LR --verbose --skinny"
    "--sort=LR --verbose --hefty"
    "--sort=SR --verbose --skinny"
    "--sort=SR --verbose --hefty"
    "--verbose"
  COMM serial mpi
  )
