
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

SET(ARGS "--filedir=${CMAKE_CURRENT_BINARY_DIR}/")

IF (SITE STREQUAL "gabriel.sandia.gov")
  SET(ARGS
    "${ARGS} --not-unit-test=CrsMatrix_int_ComplexFloat_FullMatrixComplex_UnitTest")
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_UnitTests
  SOURCES
    CrsMatrix_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi tpetratestlib
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_BlockExtraction
  SOURCES
    CrsMatrix_BlockExtraction
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  # ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  DEPLIBS tpetraext tpetratestlib
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(TpetraCrsMatrixCopyFiles
  SOURCE_FILES west0067.rua mhd1280b.cua
  EXEDEPS CrsMatrix_UnitTests
  )

ASSERT_DEFINED(Anasazi_SOURCE_DIR)
TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_CRSTests_CopyFiles1
  SOURCE_DIR ${Anasazi_SOURCE_DIR}/testmatrices
  SOURCE_FILES bcsstk14.hb bcsstk15.hb
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_InOutTests
  SOURCES
    CrsMatrix_InOutTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi tpetrainout tpetratestlib
  )

SET(TIMING_INSTALLS "")

#IF (KokkosClassic_ENABLE_TBB AND Tpetra_ENABLE_MPI)
#  TRIBITS_ADD_EXECUTABLE(
#      CRSTiming_TBB
#      SOURCES CRSTiming_TBB
#      COMM mpi
#      DEPLIBS kokkos kokkoslinalg kokkosnodeapi
#      )
#  APPEND_SET(TIMING_INSTALLS Tpetra_CRSTiming_TBB)
#ENDIF()

#IF (KokkosClassic_ENABLE_ThreadPool AND Tpetra_ENABLE_MPI)
#  TRIBITS_ADD_EXECUTABLE(
#      CRSTiming_TPI
#      SOURCES CRSTiming_TPI
#      COMM mpi
#      DEPLIBS kokkos kokkoslinalg kokkosnodeapi
#      )
#  APPEND_SET(TIMING_INSTALLS Tpetra_CRSTiming_TPI)
#ENDIF()

INSTALL(TARGETS ${TIMING_INSTALLS}
        RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}")
