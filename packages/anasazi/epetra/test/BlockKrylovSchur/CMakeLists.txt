INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Epetra_BKS_solvertest
  SOURCES cxx_main_solvertest.cpp
  ARGS 
    "--verbose"
#    "--debug"
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Epetra_BKS_test
  SOURCES cxx_main.cpp
  ARGS 
    "--verbose --exsitu"
    "--verbose --insitu"
#    "--debug"
#    "--debug --shortrun"
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Epetra_BKS_nh_test
  SOURCES cxx_main_nh.cpp
  ARGS 
    "--verbose --exsitu"
    "--verbose --insitu"
#    "--debug"
#    "--debug --shortrun"
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Epetra_BKS_QEVP_test
  SOURCES cxx_qevp.cpp
  ARGS 
    "--A-filename=linearized_qevp_A.hb --B-filename=linearized_qevp_B.hb"
  COMM serial
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestBlockKrylovSchurQEVPFiles
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/anasazi/epetra/test/BlockKrylovSchur
  SOURCE_FILES linearized_qevp_A.hb linearized_qevp_B.hb
  EXEDEPS Epetra_BKS_QEVP_test
  )

