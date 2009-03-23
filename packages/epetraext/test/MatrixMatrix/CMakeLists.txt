
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MatrixMatrix_test
  SOURCES cxx_main.cpp
  ARGS -v
  COMM serial mpi
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(MatrixMatrixCopyFiles
  DEST_FILES
    infileAB infileATBT infileAB2 infiles infileABT infileAB3 infileATB infileATB2
    C.mtx C4x4.mtx C4x12x12x4.mtx C4x6.mtx C6x4.mtx C6x6.mtx mat6x4.mtx mat6x6.mtx
    mat12x4.mtx mat4x12.mtx mat4x4.mtx mat4x6.mtx Y.mtx YTC.mtx Y_transp.mtx
    roman roman2 romancase romancase2 cdt cdt_case cdt_d.mtx cdt_m.mtx cdt_tce.mtx
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/test/MatrixMatrix
  SOURCE_PREFIX "src_"
  EXEDEPS MatrixMatrix_test
  )

  
