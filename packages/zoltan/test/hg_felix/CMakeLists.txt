####  Main Zoltan CTEST File
####  This file can be copied and modified for new Zoltan test directories.
####  It runs script zoltan/test/ctest_zoltan.pl.
####  See ctest_zoltan.pl for the list of arguments to provide to it.
####
####  To use this file in new test directories, copy it to the new test
####  directory.  In the copy, modify the number of processors to use in MPI
####  runs, ZOLTANNUMPROCS, in the SET command below.
####  Also modify the test name (e.g., substitute the new directory name
####  for "hg_felix" in all lines below.)

# Number of processors to use in MPI runs.
SET(ZOLTANNUMPROCS 5)

INCLUDE(PackageAddTest)
INCLUDE(PackageAddAdvancedTest)
INCLUDE(PackageCopyFilesToBinaryDir)

# All zdrive.inp* and hier*inp files needed for testing.
FILE(GLOB ZDRIVE_INP RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} zdrive.inp.*)
FILE(GLOB HIER_INP RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} hier*.inp)

# All accepted answer files needed for testing.
FILE(GLOB ANSWERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} answers/felix.*)


PACKAGE_ADD_ADVANCED_TEST(
  hg_felix_zoltan_parallel
  # Test Native Zoltan algorithms in parallel
  TEST_0 ${ZOLTAN_TEST_DRIVER} ${ZOLTANNUMPROCS}  Zoltan 1 ${MPI_EXEC}
  COMM mpi
  OVERALL_NUM_MPI_PROCS ${ZOLTANNUMPROCS} 
  FINAL_FAIL_REGULAR_EXPRESSION "FAILED"
  )

IF (${PACKAGE_NAME}_ENABLE_ParMETIS)
PACKAGE_ADD_ADVANCED_TEST(
  hg_felix_parmetis_parallel
  # Test ParMETIS algorithms in parallel
  TEST_0 ${ZOLTAN_TEST_DRIVER} ${ZOLTANNUMPROCS} ParMETIS 1 ${MPI_EXEC}
  COMM mpi
  OVERALL_NUM_MPI_PROCS ${ZOLTANNUMPROCS}
  FINAL_FAIL_REGULAR_EXPRESSION "FAILED"
  )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Scotch)
PACKAGE_ADD_ADVANCED_TEST(
  hg_felix_scotch_parallel
  # Test Scotch algorithms in parallel
  TEST_0 ${ZOLTAN_TEST_DRIVER} ${ZOLTANNUMPROCS} Scotch 1 ${MPI_EXEC}
  COMM mpi
  OVERALL_NUM_MPI_PROCS ${ZOLTANNUMPROCS}
  FINAL_FAIL_REGULAR_EXPRESSION "FAILED"
  )
ENDIF()

PACKAGE_ADD_ADVANCED_TEST(
  hg_felix_zoltan_serial
  # Test Native Zoltan algorithms in serial
  TEST_0 ${ZOLTAN_TEST_DRIVER} 1 Zoltan 1
  COMM serial
  OVERALL_NUM_MPI_PROCS 1
  FINAL_FAIL_REGULAR_EXPRESSION "FAILED"
  )

IF (${PACKAGE_NAME}_ENABLE_ParMETIS)
PACKAGE_ADD_ADVANCED_TEST(
  hg_felix_parmetis_serial
  # Test ParMETIS algorithms in serial
  TEST_0 ${ZOLTAN_TEST_DRIVER} 1 ParMETIS 1
  COMM serial
  OVERALL_NUM_MPI_PROCS 1
  FINAL_FAIL_REGULAR_EXPRESSION "FAILED"
  )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Scotch)
PACKAGE_ADD_ADVANCED_TEST(
  hg_felix_scotch_serial
  # Test Scotch algorithms in serial
  TEST_0 ${ZOLTAN_TEST_DRIVER} 1 Scotch 1
  COMM serial
  OVERALL_NUM_MPI_PROCS 1
  FINAL_FAIL_REGULAR_EXPRESSION "FAILED"
  )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_PaToH)
PACKAGE_ADD_ADVANCED_TEST(
  hg_felix_patoh_serial
  # Test PaToH algorithms in serial
  TEST_0 ${ZOLTAN_TEST_DRIVER} 1 PaToH 1
  COMM serial
  OVERALL_NUM_MPI_PROCS 1
  FINAL_FAIL_REGULAR_EXPRESSION "FAILED"
  )
ENDIF()


# Command that copies files to the executable directory.
PACKAGE_COPY_FILES_TO_BINARY_DIR(ZDriveHgFelixCopyFiles
  SOURCE_FILES
    felix.mtxp
    ${ZDRIVE_INP}
    ${HIER_INP}
    ${ANSWERS}
  EXEDEPS zdrive
  NOEXEPREFIX
  )
