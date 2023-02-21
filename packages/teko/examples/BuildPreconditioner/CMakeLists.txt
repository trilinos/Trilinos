ADD_SUBDIRECTORY(step1)


TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack2)
ASSERT_DEFINED(Tpetra_INST_DOUBLE)

IF (Tpetra_INST_DOUBLE)
  IF (TEKO_ENABLE_EPETRA)
    TRIBITS_ADD_EXECUTABLE(
      BuildPreconditioner_Driver_Example
      SOURCES
        example-driver.cpp
      COMM mpi
    )
  ENDIF()

  IF(${PACKAGE_NAME}_ENABLE_Ifpack2)
    TRIBITS_ADD_EXECUTABLE(
      BuildPreconditioner_Driver_Example-Belos
      SOURCES
        example-driver-belos.cpp
       COMM mpi
    )
  ENDIF()
ENDIF ()
