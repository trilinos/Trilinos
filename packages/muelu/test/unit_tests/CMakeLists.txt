TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/)
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/utils)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers/BlockedTransfers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers/BlockedSmoothers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Utils)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Smoothers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../research/regionMG/src)

#TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

#TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

SET(SOURCES "")
SET(SOURCES_INTREPID2 "")
SET(SOURCES_BLOCKED "")
SET(SOURCES_CREATEPRECONDITIONER "")
IF (${PACKAGE_NAME}_ENABLE_Experimental)
  SET(SOURCES_REGION "")
ENDIF()

APPEND_SET(SOURCES
  Aggregates.cpp
  AggregateQualityEstimateFactory.cpp
  AmalgamationFactory.cpp
  AmalgamationInfo.cpp
  BlackBoxPFactory.cpp
  CoalesceDropFactory.cpp
  CoarseMapFactory.cpp
  ClassicalPFactory.cpp
  FineLevelInputDataFactory.cpp
  GeneralGeometricPFactory.cpp
  GenericRFactory.cpp
  Hierarchy.cpp
  IndexManager.cpp
  Level.cpp
  LowPrecision.cpp
  LWGraph.cpp
  MapTransferFactory.cpp
  Memory.cpp
  MueLu_UnitTests.cpp
  MueLu_TestHelpers.cpp
  MultiVectorTransferFactory.cpp
  NullspaceFactory.cpp
  SemiCoarsenPFactoryWithSemiRestriction.cpp
  ParameterList/FactoryFactory.cpp
  ParameterList/ParameterListInterpreter.cpp
  PermutedTransferFactory.cpp
  AlgebraicPermutationStrategy.cpp
  PgPFactory.cpp
  RAPFactory.cpp
  RAPShiftFactory.cpp
  RebalanceAcFactory.cpp
  ReitzingerPFactory.cpp
  Regression.cpp
  SaPFactory.cpp
  ScaledNullspaceFactory.cpp
  SegregatedAFactory.cpp
  SemiCoarsenPFactory.cpp
  Smoothers/SmootherFactory.cpp
  StructuredAggregationFactory.cpp
  StructuredLineDetection.cpp
  TentativePFactory.cpp
  ThresholdAFilterFactory.cpp
  TransPFactory.cpp
  UncoupledAggregationFactory.cpp
  UnsmooshFactory.cpp
  UserData/CreateXpetraPreconditioner.cpp
  Utilities.cpp
  VariableContainer.cpp
  VariableDofLaplacianFactory.cpp
  MueLu_CoupledRBMFactory.cpp
  PermutationFactory.cpp
  RigidBodyModeFactory.cpp
  ML2MueLuParameterTranslator.cpp
)

APPEND_SET(SOURCES_INTREPID2
  MueLu_UnitTests.cpp
  MueLu_TestHelpers.cpp
)

APPEND_SET(SOURCES_BLOCKED
  MueLu_UnitTests.cpp
  MueLu_TestHelpers.cpp
)

IF (${PACKAGE_NAME}_ENABLE_Experimental)
  APPEND_SET(SOURCES_REGION
    MueLu_UnitTests.cpp
    MueLu_TestHelpers.cpp
  )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(UnitTestsUserData_cp
  SOURCE_DIR ${MueLu_SOURCE_DIR}/test/unit_tests/UserData
  SOURCE_FILES test.xml
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/UserData
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(UnitTestsTestMatrices_cp
  SOURCE_DIR ${MueLu_SOURCE_DIR}/test/unit_tests/TestMatrices
  SOURCE_FILES cd2dx.mat
  SOURCE_FILES cd2dy.mat
  SOURCE_FILES aniso2dx.mat
  SOURCE_FILES aniso2dy.mat
  SOURCE_FILES iso2d.mat
  SOURCE_FILES beam.mm
  SOURCE_FILES filter.mm
  SOURCE_FILES fuego0.mm
  SOURCE_FILES fuego1.mm
  SOURCE_FILES nonsym.mm
  SOURCE_FILES SaP_constrainTest_P.mat
  SOURCE_FILES semiRblkTestver.mm
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/TestMatrices
)

IF (${PACKAGE_NAME}_ENABLE_Epetra)

IF( HAVE_${PACKAGE_NAME_UC}_DEPRECATED_CODE)
APPEND_SET(SOURCES
  ParameterList/MLParameterListInterpreter.cpp
)
ENDIF()

APPEND_SET(SOURCES
  ParameterList/CreateSublists.cpp
)
ENDIF()

APPEND_SET(SOURCES_BLOCKED
  BlockedCoarseMapFactory.cpp
  BlockedPFactory.cpp
  BlockedRAPFactory.cpp
  InverseApproximationFactory.cpp
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
  ${PACKAGE_NAME}_ENABLE_Belos
  ${PACKAGE_NAME}_ENABLE_Zoltan
  )
APPEND_SET(SOURCES Adapters/TpetraOperatorAdapter.cpp)


IF(${PACKAGE_NAME}_ENABLE_Ifpack2)
  APPEND_SET(SOURCES Smoothers/Ifpack2Smoother.cpp)
  APPEND_SET(SOURCES_BLOCKED Smoothers/BlockedSmoother.cpp)
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Ifpack)
  APPEND_SET(SOURCES Smoothers/IfpackSmoother.cpp)
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Belos)
  APPEND_SET(SOURCES Smoothers/BelosSmoother.cpp)
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Stratimikos AND ${PACKAGE_NAME}_ENABLE_Thyra)
  APPEND_SET(SOURCES Smoothers/StratimikosSmoother.cpp)
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Amesos)
  APPEND_SET(SOURCES Smoothers/AmesosSmoother.cpp)
  IF(NOT ${PACKAGE_NAME}_ENABLE_Amesos2)
    APPEND_SET(SOURCES_BLOCKED Smoothers/BlockedDirectSolver.cpp)
  ENDIF()
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Amesos2)
  APPEND_SET(SOURCES Smoothers/Amesos2Smoother.cpp)
  APPEND_SET(SOURCES_BLOCKED Smoothers/BlockedDirectSolver.cpp)
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_AmgX)
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
    (${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2))
  SET(${PACKAGE_NAME}_ENABLE_CreatePreconditionerTests ON)
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_CreatePreconditionerTests OFF)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_CreatePreconditionerTests)
  APPEND_SET(SOURCES_CREATEPRECONDITIONER
    MueLu_UnitTests.cpp
    MueLu_TestHelpers.cpp
  )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_CreatePreconditionerTests)
  APPEND_SET(SOURCES_CREATEPRECONDITIONER Adapters/CreatePreconditioner.cpp)

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

IF (${PACKAGE_NAME}_ENABLE_Experimental AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)
  APPEND_SET(SOURCES_REGION
    RegionMatrix.cpp
    RegionVector.cpp
    RegionRFactory.cpp)
ENDIF()


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

  IF (${PACKAGE_NAME}_ENABLE_CreatePreconditionerTests)
    TRIBITS_ADD_EXECUTABLE(
      UnitTests_CreatePreconditioner
      SOURCES ${SOURCES_CREATEPRECONDITIONER}
      COMM serial mpi
      )
  ENDIF()

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
      UnitTests_Blocked
      NAME "UnitTestsBlockedEpetra"
      ARGS "--linAlgebra=Epetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 4
      COMM mpi
      )

    IF (${PACKAGE_NAME}_ENABLE_CreatePreconditionerTests)
      TRIBITS_ADD_TEST(
        UnitTests_CreatePreconditioner
        NAME "UnitTestsCreatePreconditionerEpetra"
        ARGS "--linAlgebra=Epetra"
        PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
        NUM_MPI_PROCS 1
        COMM serial mpi
        )

      TRIBITS_ADD_TEST(
        UnitTests_CreatePreconditioner
        NAME "UnitTestsCreatePreconditionerEpetra"
        ARGS "--linAlgebra=Epetra"
        PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
        NUM_MPI_PROCS 4
        COMM mpi
        )
    ENDIF()

  ENDIF() # Epetra / EpetraExt

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
    RUN_SERIAL
    )

  TRIBITS_ADD_TEST(
    UnitTests_Intrepid2
    NAME "UnitTestsIntrepid2Tpetra"
    ARGS "--linAlgebra=Tpetra"
    PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
    NUM_MPI_PROCS 4
    COMM mpi
    RUN_SERIAL
    )

  IF (${PACKAGE_NAME}_ENABLE_CreatePreconditionerTests)
    TRIBITS_ADD_TEST(
      UnitTests_CreatePreconditioner
      NAME "UnitTestsCreatePreconditionerTpetra"
      ARGS "--linAlgebra=Tpetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      UnitTests_CreatePreconditioner
      NAME "UnitTestsCreatePreconditionerTpetra"
      ARGS "--linAlgebra=Tpetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 4
      COMM mpi
      )
  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Experimental)
    APPEND_SET(SOURCES_REGION Galeri.cpp)

    TRIBITS_ADD_EXECUTABLE(
      UnitTests_Region
      SOURCES ${SOURCES_REGION}
      COMM serial mpi
    )

    TRIBITS_ADD_TEST(
      UnitTests_Region
      NAME "UnitTestsRegion"
      ARGS "--linAlgebra=Tpetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      UnitTests_Region
      NAME "UnitTestsRegion"
      ARGS "--linAlgebra=Tpetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 4
      COMM mpi
      )
  ENDIF()



ADD_SUBDIRECTORY(ParameterList/FactoryFactory/)
ADD_SUBDIRECTORY(ParameterList/ParameterListInterpreter/)

IF (${PACKAGE_NAME}_ENABLE_Epetra)
  ADD_SUBDIRECTORY(ParameterList/MLParameterListInterpreter/)
  ADD_SUBDIRECTORY(ParameterList/CreateSublists/)
ENDIF()
