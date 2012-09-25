
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/utils)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers/BlockedTransfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers/BlockedSmoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)

#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

#TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

SET(SOURCES "")

APPEND_SET(SOURCES
  MueLu_UnitTests.cpp
  MueLu_TestHelpers.cpp
  Level.cpp
  Hierarchy.cpp
  RAPFactory.cpp
  SaPFactory.cpp
  PgPFactory.cpp
  TransPFactory.cpp
  Needs.cpp
  TentativePFactory.cpp
  Aggregates.cpp
  CoalesceDropFactory.cpp
  Memory.cpp
  #UCAggregationFactory.cpp
  TwoKeyMap.cpp
  GenericRFactory.cpp
  ThresholdAFilterFactory.cpp
  Smoothers/SmootherFactory.cpp
  Smoothers/MueLu_TestHelpersSmoothers.cpp
  Utilities
  VariableContainer.cpp
  SubBlockAFactory.cpp
  ParameterList/FactoryFactory.cpp
  ParameterList/ParameterListInterpreter.cpp
  MultiVectorTransferFactory.cpp
)

### Tests that require other Trilinos packages

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Ifpack
  ${PACKAGE_NAME}_ENABLE_Ifpack2
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_EpetraExt
  ${PACKAGE_NAME}_ENABLE_Tpetra
  ${PACKAGE_NAME}_ENABLE_Belos
  ${PACKAGE_NAME}_ENABLE_Zoltan
  )

IF (${PACKAGE_NAME}_ENABLE_Ifpack2)
  APPEND_SET(SOURCES Smoothers/Ifpack2Smoother.cpp)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Ifpack)
  APPEND_SET(SOURCES Smoothers/IfpackSmoother.cpp)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Amesos)
  APPEND_SET(SOURCES Smoothers/AmesosSmoother.cpp)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Amesos2)
  APPEND_SET(SOURCES Smoothers/Amesos2Smoother.cpp)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Belos)
  APPEND_SET(SOURCES BelosAdapters.cpp)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Zoltan AND TPL_ENABLE_MPI)
  APPEND_SET(SOURCES Zoltan.cpp
              Repartition.cpp
              #PermutedTransferFactory.cpp  # TODO CHECK ME!
  )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Galeri) #TMP
  TRIBITS_ADD_EXECUTABLE(
    UnitTests
    SOURCES ${SOURCES}
    COMM serial mpi
    DEPLIBS muelu muelu-xpetra
    )
  
  IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt
      AND ${PACKAGE_NAME}_ENABLE_Ifpack  AND ${PACKAGE_NAME}_ENABLE_Amesos)
    
    TRIBITS_ADD_TEST(
      UnitTests
      NAME "UnitTestsEpetra"
      ARGS "--linAlgebra=0"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 1
      COMM serial mpi
      )
    
    TRIBITS_ADD_TEST(
      UnitTests
      NAME "UnitTestsEpetra"
      ARGS "--linAlgebra=0"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 4
      COMM mpi
      )
    
  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)
    
    TRIBITS_ADD_TEST(
      UnitTests
      NAME "UnitTestsTpetra"
      ARGS "--linAlgebra=1"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 1
      COMM serial mpi
      )
    
    TRIBITS_ADD_TEST(
      UnitTests
      NAME "UnitTestsTpetra"
      ARGS "--linAlgebra=1"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 4
      COMM mpi
      )
    
  ENDIF()

ENDIF() # Galeri TMP

ADD_SUBDIRECTORY(ParameterList/FactoryFactory/)
ADD_SUBDIRECTORY(ParameterList/ParameterListInterpreter/)
