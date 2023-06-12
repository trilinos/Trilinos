
IF(${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_ThyraEpetraAdapters)
TRIBITS_ADD_EXECUTABLE(
  Thyra_MVOPTester_example
  SOURCES MVOPTesterEx.cpp
  COMM serial mpi
)
ENDIF()
