
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Galeri)

TRIBITS_ADD_EXECUTABLE(
  BlockGmres_Tpetra_File_Ex
  SOURCES BlockGmresTpetraExFile.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  PseudoBlockGmres_Tpetra_File_Ex
  SOURCES PseudoBlockGmresTpetraExFile.cpp
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_CopyExampleBlockGmresFiles
  SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
  SOURCE_FILES orsirr1.hb identity.mtx
  EXEDEPS BlockGmres_Tpetra_File_Ex
  )

IF (${PACKAGE_NAME}_ENABLE_Galeri AND ${PACKAGE_NAME}_ENABLE_Xpetra)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_BlockGmres_Galeri_Ex
    SOURCES BlockGmresTpetraGaleriEx.cpp
    COMM serial mpi
    ARGS "--verbose --nx=10 --ortho=DGKS"
         "--verbose --nx=10 --ortho=ICGS"
         "--verbose --nx=10 --ortho=IMGS"
    STANDARD_PASS_OUTPUT
    )

ENDIF (${PACKAGE_NAME}_ENABLE_Galeri AND ${PACKAGE_NAME}_ENABLE_Xpetra)
