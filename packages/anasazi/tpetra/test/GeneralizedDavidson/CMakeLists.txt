
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_GeneralizedDavidson_Lap_test
  SOURCES cxx_main_lap.cpp
  ARGS
  COMM serial mpi
  )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_GeneralizedDavidson_test
    SOURCES cxx_main.cpp
    ARGS
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_GeneralizedDavidson_CopyFiles
    SOURCE_DIR ${PACKAGE_SOURCE_DIR}/testmatrices
    SOURCE_FILES mhd1280b.cua
    EXEDEPS Tpetra_GeneralizedDavidson_test
    )
ENDIF()


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_GeneralizedDavidson_nh_test
  SOURCES cxx_main_nh.cpp
  ARGS
    "--verbose"
#    "--debug"
  COMM serial mpi
  )
