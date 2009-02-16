INCLUDE(PackageAddExecutableAndTest)
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

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)

IF (${PACKAGE_NAME}_ENABLE_Amesos)
  IF (${PACKAGE_NAME}_ENABLE_EpetraExt)

  PACKAGE_ADD_EXECUTABLE(
    Epetra_BKS_QEVP_test
    SOURCES cxx_qevp.cpp
    COMM serial mpi
    )

  PACKAGE_ADD_TEST(
    Epetra_BKS_QEVP_test
    ARGS 
      "--A-filename=linearized_qevp_A.hb --B-filename=linearized_qevp_B.hb"
    COMM serial
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestBlockKrylovSchurQEVPFiles
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/anasazi/testmatrices
    SOURCE_FILES linearized_qevp_A.hb linearized_qevp_B.hb
    EXEDEPS Epetra_BKS_QEVP_test
    )
ENDIF(${PACKAGE_NAME}_ENABLE_EpetraExt)
ENDIF(${PACKAGE_NAME}_ENABLE_Amesos)
