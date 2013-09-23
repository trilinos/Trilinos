INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Epetra
  )

ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_EpetraExt
  )

ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Galeri
  )

ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  )

ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_AztecOO
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  laplace2d
  SOURCES laplace2d.cpp
  NUM_MPI_PROCS 2
  COMM mpi serial
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(laplace2d_cp
  SOURCE_FILES hands-on.sh MueLu_Agg2VTK.py xml/s2a.xml xml/s2b.xml xml/s4a.xml
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  recirc2d
  SOURCES recirc2d.cpp
  NUM_MPI_PROCS 2
  COMM mpi serial
  )
  
TRIBITS_COPY_FILES_TO_BINARY_DIR(recirc2d_cp
  SOURCE_FILES hands-on.sh xml/s3a.xml xml/s3b.xml xml/s3c.xml
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  recirc2d_api
  SOURCES recirc2d_api.cpp
  NUM_MPI_PROCS 2
  COMM mpi serial
  )  
  