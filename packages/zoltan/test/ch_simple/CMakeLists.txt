
INCLUDE(Trilinos_Add_Test)
INCLUDE(Copy_Files_To_Binary_Dir)


TRILINOS_ADD_TEST(
  ../../src/driver/zdrive
  NOEXEPREFIX
  NAME zdrive_ch_simple_rcb
  ARGS zdrive.inp.rcb
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "fatal: "
  )
# 2008/11/18: rabartl: ToDo: Above, we are *not* running the MPI test because
# it fails.  We need to talk with the Zoltan developers about this!

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ParMETIS)
IF (${PACKAGE_NAME}_ENABLE_ParMETIS)
  TRILINOS_ADD_TEST(
    ../../src/driver/zdrive
    NOEXEPREFIX
    NAME zdrive_ch_simple_partkway-v3
    ARGS zdrive.inp.partkway-v3
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "fatal: "
    )
ENDIF()


COPY_FILES_TO_BINARY_DIR(ZDriveChSimpleCopyFiles
  SOURCE_FILES
    zdrive.inp.partkway-v3
    simple.graph
    zdrive.inp.rcb
  EXEDEPS zdrive
  NOEXEPREFIX
  )
