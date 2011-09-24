INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_COPY_FILES_TO_BINARY_DIR(zoltan2_test_data_files
  SOURCE_FILES cage10.mtx  diag500_4.mtx  simple.mtx
  NOEXEPREFIX
)
