
IF(TPL_ENABLE_MPI AND ${PACKAGE_NAME}_ENABLE_Tpetra AND
   ${PACKAGE_NAME}_ENABLE_Amesos2 AND ${PACKAGE_NAME}_ENABLE_Belos AND
   ${PACKAGE_NAME}_ENABLE_Intrepid2 AND ${PACKAGE_NAME}_ENABLE_MueLu AND
   ${PACKAGE_NAME}_ENABLE_STKIO AND ${PACKAGE_NAME}_ENABLE_STKMesh AND
   NOT Kokkos_ENABLE_CUDA AND
   NOT Kokkos_ENABLE_HIP AND
   NOT Kokkos_ENABLE_SYCL)
 message(STATUS "TrilinosCouplings: Building Example_Poisson_STK")
 TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_Poisson_STK
    SOURCES example_Poisson_stk.cpp
    TrilinosCouplings_TpetraIntrepidPoissonExample_SolveWithBelos.cpp
    NUM_MPI_PROCS 1
    )
  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyMeshFilesSTK
  SOURCE_FILES unit_cube_10int_hex.exo
               unit_cube_5int_tet.exo
	       unit_cube_int_tri.exo
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Belos AND
   ${PACKAGE_NAME}_ENABLE_Intrepid2 AND
   ${PACKAGE_NAME}_ENABLE_Pamgen AND
   ${PACKAGE_NAME}_ENABLE_Tpetra AND
   NOT Kokkos_ENABLE_CUDA AND
   NOT Kokkos_ENABLE_HIP AND
   NOT Kokkos_ENABLE_SYCL)
  TRIBITS_ADD_EXECUTABLE(
    Example_Poisson2D_pn_tpetra
    SOURCES example_Poisson2D_pn_tpetra.cpp
    TrilinosCouplings_Pamgen_Utils.cpp
    TrilinosCouplings_IntrepidPoissonExampleHelpers.cpp
    COMM mpi
    )

  TRIBITS_ADD_TEST(
    Example_Poisson2D_pn_tpetra
    NAME "Poisson2D_pn=4_tpetra"
    ARGS "--mesh=Poisson2D_p4.xml --solver=muelu_ho.xml"
    COMM mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(TestCopyPoisson2DFiles
    SOURCE_FILES Poisson2D_p4.xml muelu_ho.xml
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )

  IF(${PACKAGE_NAME}_ENABLE_MueLu)
    TRIBITS_ADD_EXECUTABLE(
      Example_Maxwell_Tpetra
      SOURCES example_Maxwell_Tpetra.cpp
      TrilinosCouplings_Pamgen_Utils.cpp
      COMM serial mpi
    )

    TRIBITS_ADD_TEST(
      Example_Maxwell_Tpetra
      NAME Example_Maxwell_Tpetra_MueLu
      ARGS "--solverName=MueLu"
      NUM_MPI_PROCS 1
      COMM mpi serial
    )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(TestCopyMaxwellInputFile
      SOURCE_FILES Maxwell.xml stratimikos.xml
      DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )
  ENDIF()
ENDIF()


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
