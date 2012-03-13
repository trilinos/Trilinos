INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(TribitsLibraryMacros)

SET(TESTDATA "")
APPEND_SET(TESTDATA diag500_4.mtx)     
APPEND_SET(TESTDATA USAir97.mtx USAir97_coord.mtx)
APPEND_SET(TESTDATA commanche_dual.mtx  commanche_dual_coord.mtx)
APPEND_SET(TESTDATA simple.mtx simple_coord.mtx)

TRIBITS_COPY_FILES_TO_BINARY_DIR(zoltan2_test_data_files
  SOURCE_FILES ${TESTDATA}
  NOEXEPREFIX
)
