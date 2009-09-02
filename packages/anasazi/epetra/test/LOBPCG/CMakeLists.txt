INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Epetra_LOBPCG_solvertest
  SOURCES cxx_main_solvertest.cpp
  ARGS 
    "--verbose"
#    "--debug"
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Epetra_LOBPCG_simpletest
  SOURCES cxx_simple.cpp
  ARGS 
    "--verbose"
#    "--debug"
#    "--debug --shortrun"
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Epetra_LOBPCG_test
  SOURCES cxx_main.cpp
  ARGS 
    "--verbose"
#    "--debug"
#    "--debug --shortrun"
#    "--debug --testrecovery --nofullortho, --debug --testrecovery --fullortho"
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Epetra_LOBPCG_auxtest
  SOURCES cxx_main_auxtest.cpp
  ARGS 
    "--verbose"
#    "--debug"
  COMM serial mpi
  )
