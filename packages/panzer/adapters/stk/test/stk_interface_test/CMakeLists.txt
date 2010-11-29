INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  tSTKInterface
  SOURCES tSTKInterface.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  tSquareQuadMeshFactory
  SOURCES tSquareQuadMeshFactory.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  tSTK_IO
  SOURCES tSTK_IO.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  COMM serial mpi
  )

#PACKAGE_ADD_EXECUTABLE_AND_TEST(
#  tExodusReaderFactory
#  SOURCES tExodusReaderFactory.cpp ${UNIT_TEST_DRIVER}
#  NUM_MPI_PROCS 1
#  COMM serial mpi
#  )

#PACKAGE_COPY_FILES_TO_BINARY_DIR(stk_interface_copy
#  SOURCE_FILES basic.gen
#  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/meshes"
#  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/meshes"
#  EXEDEPS tExodusReaderFactory
#  )
