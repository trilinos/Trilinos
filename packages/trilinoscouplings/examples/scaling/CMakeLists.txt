# Just build this executable; don't run it as a test.
IF (${PACKAGE_NAME}_ENABLE_Pamgen AND
    ${PACKAGE_NAME}_ENABLE_SEACASExodus AND
    ${PACKAGE_NAME}_ENABLE_SEACASNemesis)
  message("-- TrilinosCouplings: Building rtpamgen")
   TRIBITS_ADD_EXECUTABLE(
     rtpamgen
     SOURCES rtpamgen.cpp
     COMM serial mpi
     )
   TRIBITS_COPY_FILES_TO_BINARY_DIR(TestRtpamgen
     SOURCE_FILES test.rtp
     SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
     DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
     )


  TRIBITS_ADD_TEST(
   rtpamgen
   NAME "RTpamgen_test"
   ARGS "--file=test.rtp"
   COMM serial mpi
   NUM_MPI_PROCS 1
   )


ENDIF()
