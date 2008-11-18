
INCLUDE(Trilinos_Add_Executable_And_Test)
INCLUDE(Copy_Files_To_Binary_Dir)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_UnitTests
  SOURCES
    CrsMatrix_UnitTests
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  ARGS "--filedir=${CURRENT_BINARY_DIR}"
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )


SET(INPUT_FILES west0067.rua mhd1280b.cua)
COPY_FILES_TO_BINARY_DIR(TpetraCrsMatrixCopyFiles
  SOURCE_FILES ${INPUT_FILES}
  DEST_FILES ${INPUT_FILES}
  )
