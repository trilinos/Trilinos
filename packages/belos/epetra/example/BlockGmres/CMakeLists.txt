
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
IF (${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_EpetraExt)
    
  TRIBITS_ADD_EXECUTABLE(
    BlockGmres_Epetra_File_Ex
    SOURCES BlockGmresEpetraExFile.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_EXECUTABLE(
    PseudoBlockGmres_Epetra_File_Ex
    SOURCES PseudoBlockGmresEpetraExFile.cpp
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyExampleBlockGmresFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/epetra/example/BlockGmres
    SOURCE_FILES orsirr1.hb identity.mtx
    EXEDEPS BlockGmres_Epetra_File_Ex
    )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF(${PACKAGE_NAME}_ENABLE_Ifpack)

    TRIBITS_ADD_EXECUTABLE(
      BlockPrecGmres_Epetra_File_Ex
      SOURCES BlockPrecGmresEpetraExFile.cpp
      COMM serial mpi
      )

    TRIBITS_ADD_EXECUTABLE(
      PseudoBlockPrecGmres_Epetra_File_Ex
      SOURCES PseudoBlockPrecGmresEpetraExFile.cpp
      COMM serial mpi
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      BlockFlexPrecGmres_Epetra_File_Ex
      SOURCES BlockFlexGmresEpetraExFile.cpp
      COMM serial mpi
      ARGS "--verbose --filename=orsirr1.hb"
           "--verbose --filename=identity.mtx --use-rhs=true"
      STANDARD_PASS_OUTPUT
      )

    TRIBITS_ADD_EXECUTABLE(
      BlockGmresPoly_Epetra_File_Ex
      SOURCES BlockGmresPolyEpetraExFile.cpp
      COMM serial mpi
      )

  ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_EpetraExt)
