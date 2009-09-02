INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Epetra_BlockDavidson_solvertest
  SOURCES cxx_main_solvertest.cpp
  ARGS 
    "--verbose"
#    "--debug"
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Epetra_BlockDavidson_test
  SOURCES cxx_main.cpp
  ARGS 
    "--verbose"
    "--verbose --exsitu --locking"
    "--verbose --exsitu --nolocking"
    "--verbose --insitu --locking"
    "--verbose --insitu --nolocking"
    "--verbose --insitu --rblocks=4 --locking"
    "--verbose --exsitu --rblocks=4 --locking"
    "--verbose --exsitu --rblocks=4 --nolocking"
    "--verbose --insitu --rblocks=4 --nolocking"
#    "--debug"
#    "--debug --shortrun"
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Epetra_BlockDavidson_auxtest
  SOURCES cxx_main_auxtest.cpp
  ARGS 
    "--verbose"
#    "--debug"
  COMM serial mpi
  )
