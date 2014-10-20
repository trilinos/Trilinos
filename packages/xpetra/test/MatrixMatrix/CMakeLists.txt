
IF (NOT Xpetra_Epetra_NO_32BIT_GLOBAL_INDICES)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MatrixMatrix_UnitTests
    SOURCES
      MatrixMatrix_UnitTests.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF()


IF (NOT Xpetra_Epetra_NO_64BIT_GLOBAL_INDICES)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MatrixMatrix_UnitTests_LL
    SOURCES
      MatrixMatrix_UnitTests.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    DEFINES -DXPETRA_TEST_USE_LONGLONG_GO
    )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(MatrixMatrix_UnitTests_cp
 SOURCE_FILES
   A.mat
   B.mat
   AB.mat
   AtB.mat
   ABt.mat
   AtBt.mat
)
