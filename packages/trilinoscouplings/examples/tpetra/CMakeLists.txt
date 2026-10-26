
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Tpetra)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Zoltan2Core)

IF (${PACKAGE_NAME}_ENABLE_Tpetra)

  IF (${PACKAGE_NAME}_ENABLE_Zoltan2Core)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tpetra_Zoltan2_LPTrans_Ex
      SOURCES Tpetra_Zoltan2_LPTrans_Ex.cpp
      NUM_MPI_PROCS 4
      COMM mpi
      ARGS ""
      STANDARD_PASS_OUTPUT
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tpetra_Zoltan2_LPTrans_Singleton_Ex
      SOURCES Tpetra_Zoltan2_LPTrans_Singleton_Ex.cpp
      NUM_MPI_PROCS 2
      COMM mpi
      STANDARD_PASS_OUTPUT
    )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(
      Zoltan2LPTransFileCopies
      SOURCE_FILES
        SF2_LHS_Original.mm
        SF2_LHS_Reduced.mm
        SF2_Matrix_Original.mm
        SF2_RHS_Original.mm
    )

  ENDIF()

ENDIF()
