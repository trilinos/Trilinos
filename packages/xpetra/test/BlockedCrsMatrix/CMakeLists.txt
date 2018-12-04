
#TODO: re-enable
IF (${PACKAGE_NAME}_ENABLE_Epetra)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    BlockedCrsMatrix_UnitTests
    SOURCES
      BlockedCrsMatrix_UnitTests
      ../Xpetra_UnitTests
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

IF (${PACKAGE_NAME}_ENABLE_Thyra)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ThyraBlockedOperator_UnitTests
    SOURCES
      ThyraBlockedOperator_UnitTests
      ../Xpetra_UnitTests
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF()

  TRIBITS_COPY_FILES_TO_BINARY_DIR(BlockedCrsMatrix_UnitTests_copyfiles
    SOURCE_FILES
      A.mat
      P.mat
    )
ENDIF()
