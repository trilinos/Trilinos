INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

IF (Trilinos_ENABLE_CXX11)
  PACKAGE_COPY_FILES_TO_BINARY_DIR(MultiPrecInputMatrix
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/anasazi/testmatrices
      SOURCE_FILES bcsstk14.hb
      )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(MultiPrecInputStack
      SOURCE_FILES double.xml double_float.xml qd_dd.xml qd_dd_double.xml
      )

  IF (Tpetra_ENABLE_QD) 
    PACKAGE_ADD_EXECUTABLE_AND_TEST(MultiPrecExample_qd_dd_double
      SOURCES MultiPrecExample_qd_dd_double
      ARGS "--matrix-file=bcsstk14.hb --param-file=qd_dd_double.xml"
      COMM serial mpi
      STANDARD_PASS_OUTPUT
    )
  ENDIF()

  PACKAGE_ADD_EXECUTABLE_AND_TEST(MultiPrecExample_double
    SOURCES MultiPrecExample_double
    ARGS "--matrix-file=bcsstk14.hb --param-file=double.xml"
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(MultiPrecExample_double_float
    SOURCES MultiPrecExample_double_float
    ARGS "--matrix-file=bcsstk14.hb --param-file=double_float.xml"
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

ENDIF()
