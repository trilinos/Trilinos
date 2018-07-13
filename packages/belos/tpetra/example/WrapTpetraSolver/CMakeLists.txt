
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Tpetra)
IF(${PACKAGE_NAME}_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE(
    wrapTpetraSolver
    SOURCES wrapTpetraSolver.cpp
    COMM serial mpi
    )
ENDIF()
