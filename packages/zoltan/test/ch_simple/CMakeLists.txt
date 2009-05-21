
INCLUDE(PackageAddTest)
INCLUDE(PackageAddAdvancedTest)
INCLUDE(PackageCopyFilesToBinaryDir)

# All zdrive.inp* and hier*inp files needed for testing.
FILE(GLOB ZDRIVE_INP RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} zdrive.inp.*)
FILE(GLOB HIER_INP RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} hier*.inp)

# All accepted answer files needed for testing.
FILE(GLOB ANSWERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} answers/*.*)

PACKAGE_ADD_ADVANCED_TEST(
  ch_simple_zoltan_parallel
  # Test Native Zoltan algorithms in parallel
  TEST_0 CMND ../ctest_zoltan.pl ARGS 4 Zoltan 1 ${MPI_EXEC}
  COMM mpi
  OVERALL_NUM_MPI_PROCS 4
  FINAL_FAIL_REGULAR_EXPRESSION "FAILED"
  )

IF (${PACKAGE_NAME}_ENABLE_ParMETIS)
PACKAGE_ADD_ADVANCED_TEST(
  ch_simple_parmetis_parallel
  # Test ParMETIS algorithms in parallel
  TEST_0 CMND ../ctest_zoltan.pl ARGS 4 ParMETIS 1 ${MPI_EXEC}
  COMM mpi
  OVERALL_NUM_MPI_PROCS 4
  FINAL_FAIL_REGULAR_EXPRESSION "FAILED"
  )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Scotch)
PACKAGE_ADD_ADVANCED_TEST(
  ch_simple_scotch_parallel
  # Test Scotch algorithms in parallel
  TEST_0 CMND ../ctest_zoltan.pl ARGS 4 Scotch 1 ${MPI_EXEC}
  COMM mpi
  OVERALL_NUM_MPI_PROCS 4
  FINAL_FAIL_REGULAR_EXPRESSION "FAILED"
  )
ENDIF()

PACKAGE_ADD_ADVANCED_TEST(
  ch_simple_zoltan_serial
  # Test Native Zoltan algorithms in serial
  TEST_0 CMND ../ctest_zoltan.pl ARGS 1 Zoltan 1
  COMM serial
  OVERALL_NUM_MPI_PROCS 1
  FINAL_FAIL_REGULAR_EXPRESSION "FAILED"
  )

IF (${PACKAGE_NAME}_ENABLE_ParMETIS)
PACKAGE_ADD_ADVANCED_TEST(
  ch_simple_parmetis_serial
  # Test ParMETIS algorithms in serial
  TEST_0 CMND ../ctest_zoltan.pl ARGS 1 ParMETIS 1
  COMM serial
  OVERALL_NUM_MPI_PROCS 1
  FINAL_FAIL_REGULAR_EXPRESSION "FAILED"
  )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Scotch)
PACKAGE_ADD_ADVANCED_TEST(
  ch_simple_scotch_serial
  # Test Scotch algorithms in serial
  TEST_0 CMND ../ctest_zoltan.pl ARGS 1 Scotch 1
  COMM serial
  OVERALL_NUM_MPI_PROCS 1
  FINAL_FAIL_REGULAR_EXPRESSION "FAILED"
  )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_PaToH)
PACKAGE_ADD_ADVANCED_TEST(
  ch_simple_patoh_serial
  # Test PaToH algorithms in serial
  TEST_0 CMND ../ctest_zoltan.pl ARGS 1 PaToH 1
  COMM serial
  OVERALL_NUM_MPI_PROCS 1
  FINAL_FAIL_REGULAR_EXPRESSION "FAILED"
  )
ENDIF()

# Test that Ross set up.  Its only value now is to verify when the test system
# is running.  I'll delete it later.
PACKAGE_ADD_TEST(
  ../../src/driver/zdrive
  NOEXEPREFIX
  NAME ch_simple_rcb
  ARGS zdrive.inp.rcb
  COMM serial mpi
  NUM_MPI_PROCS 4
  FAIL_REGULAR_EXPRESSION "fatal: "
  )

# Test that Ross set up.  Its only value now is to verify when the test system
# is running.  I'll delete it later.
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


# Command that copies files to the executable directory.
PACKAGE_COPY_FILES_TO_BINARY_DIR(ZDriveChSimpleCopyFiles
  SOURCE_FILES
    simple.graph
    simple.coords
    ../ctest_zoltan.pl
    ${ZDRIVE_INP}
    ${HIER_INP}
    ${ANSWERS}
  EXEDEPS zdrive
  NOEXEPREFIX
  )
