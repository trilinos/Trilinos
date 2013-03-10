INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TEUCHOS_STD_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tSTKInterface
  SOURCES tSTKInterface.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tLineMeshFactory
  SOURCES tLineMeshFactory.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tSquareQuadMeshFactory
  SOURCES tSquareQuadMeshFactory.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tSquareTriMeshFactory
  SOURCES tSquareTriMeshFactory.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  tCubeHexMeshFactory
  SOURCES tCubeHexMeshFactory.cpp ${UNIT_TEST_DRIVER}
  )
TRIBITS_ADD_TEST(
  tCubeHexMeshFactory
  NUM_MPI_PROCS 2
  CATEGORIES NIGHTLY
  )

TRIBITS_ADD_EXECUTABLE(
  tCubeTetMeshFactory
  SOURCES tCubeTetMeshFactory.cpp ${UNIT_TEST_DRIVER}
  )
TRIBITS_ADD_TEST(
  tCubeTetMeshFactory
  NUM_MPI_PROCS 2
  CATEGORIES NIGHTLY
  )

# Above test takes more than 3 minutes to run on pu241 and therefore
# should not be a *BASIC* test.  It will still run in CI and Nightly
# testing.

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tSingleBlockCubeHexMeshFactory
  SOURCES tSingleBlockCubeHexMeshFactory.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 4
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tSTK_IO
  SOURCES tSTK_IO.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tExodusReaderFactory
  SOURCES tExodusReaderFactory.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tSolutionReader
  SOURCES tSolutionReader.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  SquareMeshGen 
  SOURCES SquareMeshGen.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  ScatterSpeed_Epetra
  SOURCES   tScatterSpeed_Epetra.cpp
  )

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  tUniformRef
#  SOURCES tUniformRef.cpp ${UNIT_TEST_DRIVER}
#  NUM_MPI_PROCS 2
#  COMM serial mpi
#  )

# Panzer has no direct dependency on ThreadPool, however the search
# portion of the stk library will be disabled if ThreadPool is
# disabled.  This test will fail to compile without the missing stk
# files.  So for now we will check to see if ThreadPool is
# enabled. STK will eventually remove this dependency.
IF(STK_ENABLE_ThreadPool)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    tPointLocationSearch
    SOURCES tPointLocationSearch.cpp ${UNIT_TEST_DRIVER}
    NUM_MPI_PROCS 2
    )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(stk_interface_copy
  SOURCE_FILES basic.gen basic.gen.2.0 basic.gen.2.1
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/meshes"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/meshes"
  EXEDEPS tExodusReaderFactory
  )
