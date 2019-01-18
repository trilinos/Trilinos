INCLUDE(TribitsETISupport)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# These examples need Scalar = float and double.  If explicit template
# instantiation (ETI) is enabled, we have to make sure that those
# types are in the ETI Scalar list.
ASSERT_DEFINED (TpetraCore_ENABLE_EXPLICIT_INSTANTIATION)
IF (TpetraCore_ENABLE_EXPLICIT_INSTANTIATION)
  LIST (FIND TpetraCore_ETI_SCALARS "float" floatFindIndex)
  LIST (FIND TpetraCore_ETI_SCALARS "double" doubleFindIndex)
  IF (floatFindIndex EQUAL -1 OR doubleFindIndex EQUAL -1) # -1 means not in list
    SET (floatAndDoubleEnabled OFF)
  ELSE ()
    SET (floatAndDoubleEnabled ON)
  ENDIF ()
ELSE ()
  SET (floatAndDoubleEnabled ON)
ENDIF ()

ASSERT_DEFINED(Anasazi_SOURCE_DIR)
IF (${PROJECT_NAME}_ENABLE_CXX11 AND floatAndDoubleEnabled)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(IRTR_double
    SOURCES IRTR_double
    ARGS "--matrix-file=bcsstk14.hb --param-file=irtr_double.xml"
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
      SOURCE_FILES irtr_double.xml irtr_qd.xml 
                   run_all_irtr.csh
      EXEDEPS IRTR_double
      )

  # This enable requires the QD TPL for Scalar = qd_real.  qd_real
  # must also be in the list of ETI Scalar types, if ETI is enabled.
  IF (Tpetra_ENABLE_QD AND (NOT TpetraCore_ENABLE_EXPLICIT_INSTANTIATION OR Tpetra_INST_QD_REAL))

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

  ENDIF() # Tpetra_ENABLE_QD
ENDIF()
