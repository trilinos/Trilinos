
IF(${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters)
TRIBITS_ADD_EXECUTABLE(
  Thyra_MVOPTester_example
  SOURCES MVOPTesterEx.cpp
  COMM serial mpi
)
ENDIF()
