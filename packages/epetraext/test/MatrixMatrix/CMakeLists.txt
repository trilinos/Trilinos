
INCLUDE(Trilinos_Add_Executable_And_Test)
INCLUDE(Copy_Files_To_Binary_Dir)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  MatrixMatrix_test
  SOURCES cxx_main.cpp
  COMM serial mpi
  )

COPY_FILES_TO_BINARY_DIR(EpetraExtMatrixMatrixCopyFiles
  DEST_FILES
    infileAB infileATBT infileAB2 infiles infileABT infileAB3 infileATB infileATB2
    C.mtx C4x4.mtx C4x12x12x4.mtx C4x6.mtx C6x4.mtx C6x6.mtx mat6x4.mtx mat6x6.mtx
    mat12x4.mtx mat4x12.mtx mat4x4.mtx mat4x6.mtx Y.mtx YTC.mtx Y_transp.mtx
    roman roman2 romancase romancase2 cdt cdt_case cdt_d.mtx cdt_m.mtx cdt_tce.mtx
  SOURCE_DIR ${${PACKAGE_NAME}_SOURCE_DIR}/test/MatrixMatrix
  SOURCE_PREFIX "src_"
  )

  
