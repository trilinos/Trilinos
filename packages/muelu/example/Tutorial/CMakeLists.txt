INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Epetra
  )

ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_EpetraExt
  )

IF(${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tutorial1
  SOURCES tutorial1.cpp
  NUM_MPI_PROCS 2
  COMM mpi serial
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(tutorial1_cp
  SOURCE_FILES tutorial1a.xml tutorial1b.xml tutorial1c.xml tutorial1d.xml tutorial1e.xml
  )

ADD_SUBDIRECTORY(hands-on)

ENDIF()