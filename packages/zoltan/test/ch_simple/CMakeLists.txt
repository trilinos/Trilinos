
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

# Copy all zdrive.inp* files to executable directory.
FILE(GLOB ZDRIVE_INP RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} zdrive.inp.*)
FILE(GLOB ANSWERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} answers/*)

PACKAGE_ADD_ADVANCED_TEST(
  ch_simple_parallel
  TEST_0 CMND ../ctest_zoltan.pl ARGS 4
  COMM mpi
  OVERALL_NUM_MPI_PROCS 4
  FINAL_FAIL_REGULAR_EXPRESSION "FAILED"
  )

PACKAGE_ADD_ADVANCED_TEST(
  ch_simple_serial
  TEST_0 CMND ../ctest_zoltan.pl ARGS 1
  COMM serial
  OVERALL_NUM_MPI_PROCS 1
  FINAL_FAIL_REGULAR_EXPRESSION "FAILED"
  )

PACKAGE_ADD_TEST(
  ../../src/driver/zdrive
  NOEXEPREFIX
  NAME ch_simple_rcb
  ARGS zdrive.inp.rcb
  COMM serial mpi
  NUM_MPI_PROCS 4
  FAIL_REGULAR_EXPRESSION "fatal: "
  )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ParMETIS)
IF (${PACKAGE_NAME}_ENABLE_ParMETIS)
  PACKAGE_ADD_TEST(
    ../../src/driver/zdrive
    NOEXEPREFIX
    NAME ch_simple_partkway-v3
    ARGS zdrive.inp.partkway-v3
    COMM serial mpi
    NUM_MPI_PROCS 4
    FAIL_REGULAR_EXPRESSION "fatal: "
    )
ENDIF()


PACKAGE_COPY_FILES_TO_BINARY_DIR(ZDriveChSimpleCopyFiles
  SOURCE_FILES
    simple.graph
    simple.coords
    ../ctest_zoltan.pl
    ${ZDRIVE_INP}
    ${ANSWERS}
  EXEDEPS zdrive
  NOEXEPREFIX
  )
