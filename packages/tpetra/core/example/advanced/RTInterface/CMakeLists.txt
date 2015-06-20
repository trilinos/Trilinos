INCLUDE(TribitsETISupport)

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

IF (${PROJECT_NAME}_ENABLE_CXX11 AND floatAndDoubleEnabled)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    RTIExample
    SOURCES RTIExample
    ARGS 
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    RTIOperatorExample
    SOURCES RTIOperatorExample
    ARGS 
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    RTITutorialOperator
    SOURCES RTITutorialOperator
    ARGS 
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    RTITutorialSimple
    SOURCES RTITutorialSimple
    ARGS 
    COMM mpi
    STANDARD_PASS_OUTPUT
    )
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    RTIInlineCG
    SOURCES RTIInlineCG
    ARGS 
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
  
  # This enable requires the QD TPL for Scalar = qd_real.  qd_real
  # must also be in the list of ETI Scalar types, if ETI is enabled.
  IF (Tpetra_ENABLE_QD AND (NOT TpetraCore_ENABLE_EXPLICIT_INSTANTIATION OR Tpetra_INST_QD_REAL))
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      RTIInlineCG_qd
      SOURCES RTIInlineCG_qd
      ARGS 
      COMM serial mpi
      STANDARD_PASS_OUTPUT
      )
  ENDIF()
ENDIF()
