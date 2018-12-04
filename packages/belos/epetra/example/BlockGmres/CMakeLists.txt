
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Galeri)
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

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      BlockGmresPoly_Epetra_File_Ex
      SOURCES BlockGmresPolyEpetraExFile.cpp
      COMM serial mpi
      ARGS "--verbose --filename=orsirr1.hb --precond=none --max-degree=5 --poly-type=Arnoldi"
           "--verbose --filename=orsirr1.hb --precond=none --max-degree=5 --poly-type=Gmres"
      )

  ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_EpetraExt)

IF (${PACKAGE_NAME}_ENABLE_Galeri)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    BlockGmres_Galeri_Ex
    SOURCES BlockGmresEpetraGaleriEx.cpp
    COMM serial mpi
    ARGS "--verbose --nx=10 --ortho=DGKS"
         "--verbose --nx=10 --ortho=ICGS"
         "--verbose --nx=10 --ortho=IMGS"
    STANDARD_PASS_OUTPUT
    )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF(${PACKAGE_NAME}_ENABLE_Ifpack)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      BlockPrecGmres_Galeri_Ex
      SOURCES BlockPrecGmresEpetraGaleriEx.cpp
      COMM serial mpi
      ARGS "--verbose --nx=10 --ortho=DGKS"
           "--verbose --nx=10 --ortho=ICGS"
           "--verbose --nx=10 --ortho=IMGS"
      STANDARD_PASS_OUTPUT
      )

  ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)

ENDIF (${PACKAGE_NAME}_ENABLE_Galeri)

