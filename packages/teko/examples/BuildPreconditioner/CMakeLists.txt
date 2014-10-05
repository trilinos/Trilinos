ADD_SUBDIRECTORY(step1)

INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
   BuildPreconditioner_Driver_Example
   SOURCES 
      example-driver.cpp
   COMM mpi
   )

IF(${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE(
     BuildPreconditioner_Driver_Example-Belos
     SOURCES 
        example-driver-belos.cpp
     COMM mpi
     )
ENDIF()
