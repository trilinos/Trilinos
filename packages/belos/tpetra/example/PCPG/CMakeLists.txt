

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Tpetra)

IF (${PACKAGE_NAME}_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE(
    PCPG_Tpetra_File_Ex
    SOURCES PCPGTpetraExFile.cpp
    COMM serial mpi
    # ARGS --verbose
    )

  TRIBITS_ADD_EXECUTABLE(
    PCPG_DenseKokkos_Tpetra_File_Ex
    SOURCES PCPGTpetraDenseKokkosExFile.cpp
    COMM serial mpi
    # ARGS --verbose
    )

ENDIF()
