INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

IF (Trilinos_ENABLE_CXX11)
  TRIBITS_COPY_FILES_TO_BINARY_DIR(MultiPrecInputMatrix
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/anasazi/testmatrices
      SOURCE_FILES bcsstk14.hb
      )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(MultiPrecInputMatrix2
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/tpetra/test/BlockCG
      SOURCE_FILES bcsstk17.rsa
      )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(MultiPrecInputStack
      SOURCE_FILES double.xml double_float.xml double_double.xml qd.xml qd_dd.xml qd_dd_double.xml
                   threaded.xml mpionly.xml
      )

  IF (Tpetra_ENABLE_QD) 
    TRIBITS_ADD_EXECUTABLE_AND_TEST(MultiPrecExample_qd
      SOURCES MultiPrecExample_qd
      ARGS "--matrix-file=bcsstk14.hb --param-file=qd.xml"
      COMM serial mpi
      STANDARD_PASS_OUTPUT
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(MultiPrecExample_qd_dd
      SOURCES MultiPrecExample_qd_dd
      ARGS "--matrix-file=bcsstk14.hb --param-file=qd_dd.xml"
      COMM serial mpi
      STANDARD_PASS_OUTPUT
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(MultiPrecExample_qd_dd_double
      SOURCES MultiPrecExample_qd_dd_double
      ARGS "--matrix-file=bcsstk14.hb --param-file=qd_dd_double.xml"
      COMM serial mpi
      STANDARD_PASS_OUTPUT
    )
  ENDIF()

  TRIBITS_ADD_EXECUTABLE_AND_TEST(MultiPrecExample_double_double
    SOURCES MultiPrecExample_double_double
    ARGS "--matrix-file=bcsstk14.hb --param-file=double_double.xml"
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(MultiPrecExample_double
    SOURCES MultiPrecExample_double
    ARGS "--matrix-file=bcsstk14.hb --param-file=double.xml"
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(MultiPrecExample_double_float
    SOURCES MultiPrecExample_double_float
    ARGS "--matrix-file=bcsstk14.hb --param-file=double_float.xml"
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

ENDIF()
