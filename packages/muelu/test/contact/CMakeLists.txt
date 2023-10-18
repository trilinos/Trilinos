INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Tpetra)

IF (${PACKAGE_NAME}_ENABLE_Tpetra)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Contact
    SOURCE_FILES
      segregatedAgg_block_matrix.mm
      segregatedAgg_filtered_matrix.mm
      segregatedAgg_matrix.mm
      segregatedAgg_nullspace.mm
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

  # The Contact only works with real Scalar types,
  # so choose an instantiation with real scalars (if any is available)
  SET(Contact_Instantiation "")

  IF(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
    SET(Contact_instantiation "DOUBLE_INT_LONGLONGINT")
  ELSEIF(HAVE_MUELU_INST_DOUBLE_INT_LONGINT)
    SET(Contact_instantiation "DOUBLE_INT_LONGINT")
  ELSEIF(HAVE_MUELU_INST_DOUBLE_INT_INT)
    SET(Contact_instantiation "DOUBLE_INT_INT")
  ELSEIF(HAVE_MUELU_INST_FLOAT_INT_INT)
    SET(Contact_instantiation "FLOAT_INT_INT")
  ENDIF()

  IF(Contact_instantiation)
    TRIBITS_ADD_TEST(
        SegregatedAggFromFile
        NAME "UnsegregatedAggFromFile"
        ARGS "--probName=segregatedAgg --nDof=54 --nDofPerNode=2 --instantiation=${Contact_instantiation}"
        COMM serial mpi
        NUM_MPI_PROCS 1
        WILL_FAIL
    )
    TRIBITS_ADD_TEST(
        SegregatedAggFromFile
        NAME "SegregatedAggFromFile"
        ARGS "--probName=segregatedAgg_filtered --nDof=54 --nDofPerNode=2 --instantiation=${Contact_instantiation}"
        COMM serial mpi
        NUM_MPI_PROCS 1
    )
    TRIBITS_ADD_TEST(
        SegregatedAggFromFactory
        NAME "SegregatedAggFromFactoryInterfaceMaps"
        ARGS "--probName=segregatedAgg --nDof=54 --nDofPerNode=2 --instantiation=${Contact_instantiation}"
        COMM serial mpi
        NUM_MPI_PROCS 1
    )
    TRIBITS_ADD_TEST(
        SegregatedAggFromFactory
        NAME "SegregatedAggFromFactoryInterfaceMaps"
        ARGS "--probName=segregatedAgg --nDof=54 --nDofPerNode=2 --instantiation=${Contact_instantiation}"
        COMM serial mpi
        NUM_MPI_PROCS 3
    )
    TRIBITS_ADD_TEST(
        SegregatedAggFromFactoryBlockedA
        NAME "SegregatedAggFromFactoryBlockedAInterfaceMaps"
        ARGS "--probName=segregatedAgg --nPrimalDof=54 --nDualDof=6 --numPrimalDofPerNode=2 --numDualDofPerNode=2 --instantiation=${Contact_instantiation}"
        COMM serial mpi
        NUM_MPI_PROCS 1
    )
    TRIBITS_ADD_TEST(
        SegregatedAggFromFactoryBlockedA
        NAME "SegregatedAggFromFactoryBlockedAInterfaceMaps"
        ARGS "--probName=segregatedAgg --nPrimalDof=54 --nDualDof=6 --numPrimalDofPerNode=2 --numDualDofPerNode=2 --instantiation=${Contact_instantiation}"
        COMM serial mpi
        NUM_MPI_PROCS 3
    )
  ENDIF()

ENDIF()