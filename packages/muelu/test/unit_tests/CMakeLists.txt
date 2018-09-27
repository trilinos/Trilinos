INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/)
INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/utils)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers/BlockedTransfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers/BlockedSmoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Utils)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Smoothers)

#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

#TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

SET(SOURCES "")
SET(SOURCES_INTREPID2 "")
SET(SOURCES_BLOCKED "")

APPEND_SET(SOURCES
  Aggregates.cpp
  BlackBoxPFactory.cpp
  CoalesceDropFactory.cpp
  CoarseMapFactory.cpp
# CoupledAggregationFactory.cpp
  FineLevelInputDataFactory.cpp
  GeneralGeometricPFactory.cpp
  GenericRFactory.cpp
  Hierarchy.cpp
  IndexManager.cpp
  Level.cpp
  Memory.cpp
  MueLu_UnitTests.cpp
  MueLu_TestHelpers.cpp
  MultiVectorTransferFactory.cpp
  ParameterList/FactoryFactory.cpp
  ParameterList/ParameterListInterpreter.cpp
  PermutedTransferFactory.cpp
  PgPFactory.cpp
  RAPFactory.cpp
  RAPShiftFactory.cpp
  SaPFactory.cpp
  Smoothers/SmootherFactory.cpp
  StructuredAggregationFactory.cpp
  StructuredLineDetection.cpp
  TentativePFactory.cpp
  ThresholdAFilterFactory.cpp
  TransPFactory.cpp
  UnsmooshFactory.cpp
  UserData/CreateXpetraPreconditioner.cpp
  Utilities
  VariableContainer.cpp
  VariableDofLaplacianFactory.cpp
)

APPEND_SET(SOURCES_INTREPID2
  MueLu_UnitTests.cpp
  MueLu_TestHelpers.cpp
)

APPEND_SET(SOURCES_BLOCKED
  MueLu_UnitTests.cpp
  MueLu_TestHelpers.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(UnitTestsUserData_cp
  SOURCE_DIR ${MueLu_SOURCE_DIR}/test/unit_tests/UserData
  SOURCE_FILES test.xml
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/UserData
)

IF (${PACKAGE_NAME}_ENABLE_Epetra)
APPEND_SET(SOURCES
  ParameterList/MLParameterListInterpreter.cpp
  ParameterList/CreateSublists.cpp
)
ENDIF()

APPEND_SET(SOURCES_BLOCKED
  BlockedCoarseMapFactory.cpp
  BlockedPFactory.cpp
  SchurComplementFactory.cpp
  SubBlockAFactory.cpp
)
IF (TPL_ENABLE_MPI)
  APPEND_SET(SOURCES_BLOCKED
  BlockedRepartition.cpp
  )
ENDIF()

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
IF (${PACKAGE_NAME}_ENABLE_Tpetra)
  APPEND_SET(SOURCES Adapters/TpetraOperatorAdapter.cpp)
ENDIF()


IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2)
  APPEND_SET(SOURCES Smoothers/Ifpack2Smoother.cpp)
  APPEND_SET(SOURCES_BLOCKED Smoothers/BlockedSmoother.cpp)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Ifpack)
  APPEND_SET(SOURCES Smoothers/IfpackSmoother.cpp)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Amesos)
  APPEND_SET(SOURCES Smoothers/AmesosSmoother.cpp)
  IF (NOT ${PACKAGE_NAME}_ENABLE_Tpetra OR NOT ${PACKAGE_NAME}_ENABLE_Amesos2)
    APPEND_SET(SOURCES_BLOCKED Smoothers/BlockedDirectSolver.cpp)
  ENDIF()
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Amesos2)
  APPEND_SET(SOURCES Smoothers/Amesos2Smoother.cpp)
  APPEND_SET(SOURCES_BLOCKED Smoothers/BlockedDirectSolver.cpp)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_AmgX)
  APPEND_SET(SOURCES Adapters/AmgxOperatorAdapter.cpp)
  TRIBITS_COPY_FILES_TO_BINARY_DIR(UnitTestsAmgxAdapter_cp
    SOURCE_DIR ${MueLu_SOURCE_DIR}/test/unit_tests/Adapters
    SOURCE_FILES test.json
    )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Belos)
  APPEND_SET(SOURCES Adapters/BelosAdapters.cpp)
ENDIF()

IF ((${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND ${PACKAGE_NAME}_ENABLE_Ifpack AND ${PACKAGE_NAME}_ENABLE_Amesos) OR
    (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2))
  APPEND_SET(SOURCES Adapters/CreatePreconditioner.cpp)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(UnitTestsCreatePreconditioner_cp
    SOURCE_DIR ${MueLu_SOURCE_DIR}/test/unit_tests/Adapters
    SOURCE_FILES test.xml testReuse.xml testWithRebalance.xml testPDE.xml testPDE1.xml
    )
ENDIF()



IF (${PACKAGE_NAME}_ENABLE_Intrepid2)
APPEND_SET(SOURCES_INTREPID2
  IntrepidPCoarsenFactory.cpp
)
ENDIF()


IF (${PACKAGE_NAME}_ENABLE_Zoltan AND TPL_ENABLE_MPI)
  APPEND_SET(SOURCES Zoltan.cpp
              Repartition.cpp
              #PermutedTransferFactory.cpp  # TODO CHECK ME!
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(MueLu_Repartition_cp
    SOURCE_DIR ${MueLu_SOURCE_DIR}/test/unit_tests/
    SOURCE_FILES testCoordinates.xml
  )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Galeri) #TMP
  APPEND_SET(SOURCES Galeri.cpp)
  APPEND_SET(SOURCES_INTREPID2 Galeri.cpp)
  APPEND_SET(SOURCES_BLOCKED Galeri.cpp)

  TRIBITS_ADD_EXECUTABLE(
    UnitTests
    SOURCES ${SOURCES}
    COMM serial mpi
    )

  TRIBITS_ADD_EXECUTABLE(
    UnitTests_Intrepid2
    SOURCES ${SOURCES_INTREPID2}
    COMM serial mpi
    )

  TRIBITS_ADD_EXECUTABLE(
    UnitTests_Blocked
    SOURCES ${SOURCES_BLOCKED}
    COMM serial mpi
    )

  IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt)

    TRIBITS_ADD_TEST(
      UnitTests
      NAME "UnitTestsEpetra"
      ARGS "--linAlgebra=Epetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      UnitTests
      NAME "UnitTestsEpetra"
      ARGS "--linAlgebra=Epetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 4
      COMM mpi
      )

    TRIBITS_ADD_TEST(
      UnitTests_Blocked
      NAME "UnitTestsBlockedEpetra"
      ARGS "--linAlgebra=Epetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      UnitTests
      NAME "UnitTestsBlockedEpetra"
      ARGS "--linAlgebra=Epetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 4
      COMM mpi
      )

  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Tpetra)

    TRIBITS_ADD_TEST(
      UnitTests
      NAME "UnitTestsTpetra"
      ARGS "--linAlgebra=Tpetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      UnitTests
      NAME "UnitTestsTpetra"
      ARGS "--linAlgebra=Tpetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 4
      COMM mpi
      )

    TRIBITS_ADD_TEST(
      UnitTests_Blocked
      NAME "UnitTestsBlockedTpetra"
      ARGS "--linAlgebra=Tpetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      UnitTests_Blocked
      NAME "UnitTestsBlockedTpetra"
      ARGS "--linAlgebra=Tpetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 4
      COMM mpi
      )

    TRIBITS_ADD_TEST(
      UnitTests_Intrepid2
      NAME "UnitTestsIntrepid2Tpetra"
      ARGS "--linAlgebra=Tpetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      UnitTests_Intrepid2
      NAME "UnitTestsIntrepid2Tpetra"
      ARGS "--linAlgebra=Tpetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 4
      COMM mpi
      )

  ENDIF()

ENDIF() # Galeri TMP

ADD_SUBDIRECTORY(ParameterList/FactoryFactory/)
IF (${PACKAGE_NAME}_ENABLE_Epetra)
  ADD_SUBDIRECTORY(ParameterList/MLParameterListInterpreter/)
  ADD_SUBDIRECTORY(ParameterList/ParameterListInterpreter/)
  ADD_SUBDIRECTORY(ParameterList/CreateSublists/)
ENDIF()
