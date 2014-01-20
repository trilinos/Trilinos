
INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_UnitTests
  SOURCES
    Map_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Map_UnitTests2
  SOURCES
    Map_UnitTests2
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE(
  StridedMap_UnitTests
  SOURCES
    StridedMap_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE(
  StridedMapFactory_UnitTests
  SOURCES
    StridedMapFactory_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
)

IF (${PACKAGE_NAME}_ENABLE_Epetra)

  TRIBITS_ADD_TEST(
    StridedMap_UnitTests
    NAME "StridedMapUnitTestsEpetra"
    ARGS "--linAlgebra=Epetra"
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_TEST(
    StridedMapFactory_UnitTests
    NAME "StridedMapFactoryUnitTestsEpetra"
    ARGS "--linAlgebra=Epetra"
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Tpetra)

  TRIBITS_ADD_TEST(
    StridedMap_UnitTests
    NAME "StridedMapUnitTestsTpetra"
    ARGS "--linAlgebra=Tpetra"
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_TEST(
    StridedMapFactory_UnitTests
    NAME "StridedMapFactoryUnitTestsTpetra"
    ARGS "--linAlgebra=Tpetra"
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

ENDIF()
