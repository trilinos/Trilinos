# -*- cmake -*-

IF (${PACKAGE_NAME}_ENABLE_Tpetra)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../../tpetra/src )

# Create the ValidParameters executable
TRIBITS_ADD_EXECUTABLE(
  ValidParameters
  SOURCES
    createValidParameterList.cpp
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  ValidParameters
  NAME ValidParameters
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS "--add-xsl-header"
  )

ENDIF()
