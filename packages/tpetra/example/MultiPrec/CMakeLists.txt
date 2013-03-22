INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(TribitsETISupport)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

JOIN(allnodes   "|" FALSE ${Tpetra_ETI_NODES}  )

ASSERT_DEFINED(Anasazi_SOURCE_DIR)
IF (${PROJECT_NAME}_ENABLE_CXX11)

  TRIBITS_ETI_TYPE_EXPANSION(dfneeds "S=double|float" "LO=int" "GO=int" "N=${allnodes}")
  TRIBITS_ETI_TYPE_EXPANSION(dfneeds "SOUT=float"     "SIN=double"  "LO=int" "GO=int" "N=${allnodes}")
  TRIBITS_ADD_ETI_INSTANTIATIONS(Tpetra ${dfneeds})

  TRIBITS_ADD_EXECUTABLE_AND_TEST(IRTR_double
    SOURCES IRTR_double
    ARGS "--matrix-file=bcsstk14.hb --param-file=irtr_double.xml"
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

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

  TRIBITS_COPY_FILES_TO_BINARY_DIR(MultiPrecInputMatrix
      SOURCE_DIR ${Anasazi_SOURCE_DIR}/testmatrices
      SOURCE_FILES bcsstk14.hb
      EXEDEPS IRTR_double
      )
  TRIBITS_COPY_FILES_TO_BINARY_DIR(MultiPrecInputMatrix2
      SOURCE_DIR ${Belos_SOURCE_DIR}/tpetra/test/BlockCG
      SOURCE_FILES bcsstk17.rsa
      EXEDEPS IRTR_double
      )
  TRIBITS_COPY_FILES_TO_BINARY_DIR(MultiPrecMachineFiles
      SOURCE_DIR   ${Tpetra_MACHINE_XML_FILE_DIR} 
      SOURCE_FILES threaded.xml serial.xml mpionly.xml 
      EXEDEPS IRTR_double
      )
  TRIBITS_COPY_FILES_TO_BINARY_DIR(MultiPrecInputStack
      SOURCE_FILES double.xml double_float.xml double_double.xml qd.xml qd_dd.xml qd_dd_double.xml
                   irtr_double.xml irtr_qd.xml 
                   run_all_irtr.csh
      EXEDEPS IRTR_double
      )

  IF (Tpetra_ENABLE_QD) 

    TRIBITS_ETI_TYPE_EXPANSION(qdneeds "S=qd_real|dd_real" "LO=int" "GO=int" "N=${allnodes}")
    TRIBITS_ETI_TYPE_EXPANSION(qdneeds "SOUT=double|dd_real" "SIN=qd_real" "LO=int" "GO=int" "N=${allnodes}")
    TRIBITS_ADD_ETI_INSTANTIATIONS(Tpetra ${qdneeds})

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

    TRIBITS_ADD_EXECUTABLE_AND_TEST(IRTR_qd
      SOURCES IRTR_qd
      ARGS "--matrix-file=bcsstk14.hb --param-file=irtr_qd.xml"
      COMM serial mpi
      STANDARD_PASS_OUTPUT
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(IRTR_qd_double
      SOURCES IRTR_qd_double
      ARGS "--matrix-file=bcsstk14.hb --param-file=irtr_qd.xml"
      COMM serial mpi
      STANDARD_PASS_OUTPUT
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(MultiPrecExample_qd_dd_double
      SOURCES MultiPrecExample_qd_dd_double
      ARGS "--matrix-file=bcsstk14.hb --param-file=qd_dd_double.xml"
      COMM serial mpi
      STANDARD_PASS_OUTPUT
    )

  ENDIF() # Tpetra_ENABLE_QD

ENDIF()
