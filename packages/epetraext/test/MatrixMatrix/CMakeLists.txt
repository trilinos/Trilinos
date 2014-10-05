


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMatrix_test
  SOURCES cxx_main.cpp
  ARGS -v
  COMM serial mpi
  )


TRIBITS_COPY_FILES_TO_BINARY_DIR(MatrixMatrixCopyFiles
  DEST_FILES
    infileAB infileATBT infileAB2 infiles infileABT infileAB3 infileAB4 infileATB infileATB2
    infile6x4TI6 infileI4x6x4T infilePTxI infileIxPT
    C.mtx C4x4.mtx C4x12x12x4.mtx C4x6.mtx C6x4.mtx C6x6.mtx mat6x4.mtx mat6x6.mtx
    C6x4TxI6.mtx CI4x6x4T.mtx I3.mtx I4.mtx I6.mtx I9.mtx P.mtx PT.mtx
    mat12x4.mtx mat4x12.mtx mat4x4.mtx mat4x6.mtx Y.mtx YTC.mtx Y_transp.mtx
    roman roman2 romancase romancase2 cdt cdt_case cdt_d.mtx cdt_m.mtx cdt_tce.mtx
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/MatrixMatrix
  SOURCE_PREFIX "src_"
  EXEDEPS MatrixMatrix_test
  )
