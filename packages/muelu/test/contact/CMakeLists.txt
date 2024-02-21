TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Tpetra)

IF (${PACKAGE_NAME}_ENABLE_Tpetra)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Contact
    SOURCE_FILES
      segregatedAgg_block_matrix.mm
      segregatedAgg_filtered_matrix.mm
      segregatedAgg_matrix.mm
      segregatedAgg_nullspace.mm
      segregatedAgg_interface_dof_map_MPI1.mm
      segregatedAgg_interface_dof_map_MPI3.mm
  )

  TRIBITS_ADD_EXECUTABLE(
      SegregatedAggFromFile
      SOURCES
      SegregatedAggFromFile.cpp
      contact_Helpers.cpp
      COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
      SegregatedAggFromFactory
      SOURCES
      SegregatedAggFromFactory.cpp
      contact_Helpers.cpp
      COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
      SegregatedAggFromFactoryBlockedA
      SOURCES
      SegregatedAggFromFactoryBlockedA.cpp
      contact_Helpers.cpp
      COMM serial mpi
  )

  TRIBITS_ADD_TEST(
      SegregatedAggFromFile
      NAME "UnsegregatedAggFromFile"
      ARGS "--probName=segregatedAgg --nDof=54 --nDofPerNode=2"
      COMM serial mpi
      NUM_MPI_PROCS 1
      WILL_FAIL
  )
  TRIBITS_ADD_TEST(
      SegregatedAggFromFile
      NAME "SegregatedAggFromFile"
      ARGS "--probName=segregatedAgg_filtered --nDof=54 --nDofPerNode=2"
      COMM serial mpi
      NUM_MPI_PROCS 1
  )
  TRIBITS_ADD_TEST(
      SegregatedAggFromFactory
      NAME "SegregatedAggFromFactoryInterfaceMaps"
      ARGS "--probName=segregatedAgg --nDof=54 --nDofPerNode=2"
      COMM serial mpi
      NUM_MPI_PROCS 1
  )
  TRIBITS_ADD_TEST(
      SegregatedAggFromFactory
      NAME "SegregatedAggFromFactoryInterfaceMaps"
      ARGS "--probName=segregatedAgg --nDof=54 --nDofPerNode=2"
      COMM serial mpi
      NUM_MPI_PROCS 3
  )
  TRIBITS_ADD_TEST(
      SegregatedAggFromFactoryBlockedA
      NAME "SegregatedAggFromFactoryBlockedAInterfaceMaps"
      ARGS "--probName=segregatedAgg --nPrimalDof=54 --nDualDof=6 --numPrimalDofPerNode=2 --numDualDofPerNode=2"
      COMM serial mpi
      NUM_MPI_PROCS 1
  )
  TRIBITS_ADD_TEST(
      SegregatedAggFromFactoryBlockedA
      NAME "SegregatedAggFromFactoryBlockedAInterfaceMaps"
      ARGS "--probName=segregatedAgg --nPrimalDof=54 --nDualDof=6 --numPrimalDofPerNode=2 --numDualDofPerNode=2"
      COMM serial mpi
      NUM_MPI_PROCS 3
  )

ENDIF()