ADD_SUBDIRECTORY(step1)


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack2)
ASSERT_DEFINED(Tpetra_INST_DOUBLE)

IF (Tpetra_INST_DOUBLE)
  TRIBITS_ADD_EXECUTABLE(
    BuildPreconditioner_Driver_Example
    SOURCES 
      example-driver.cpp
     COMM mpi
  )

  IF(${PACKAGE_NAME}_ENABLE_Ifpack2)
    TRIBITS_ADD_EXECUTABLE(
      BuildPreconditioner_Driver_Example-Belos
      SOURCES 
        example-driver-belos.cpp
       COMM mpi
    )
  ENDIF()
ENDIF ()
