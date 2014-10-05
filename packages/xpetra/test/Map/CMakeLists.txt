
INCLUDE(TribitsAddExecutableAndTest)

IF (NOT Xpetra_Epetra_NO_32BIT_GLOBAL_INDICES)

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
ENDIF()

IF (NOT Xpetra_Epetra_NO_64BIT_GLOBAL_INDICES)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Map_UnitTests_LL
    SOURCES
      Map_UnitTests
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    DEFINES -DXPETRA_TEST_USE_LONGLONG_GO
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Map_UnitTests2_LL
    SOURCES
      Map_UnitTests2
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    DEFINES -DXPETRA_TEST_USE_LONGLONG_GO
    )

  TRIBITS_ADD_EXECUTABLE(
    StridedMap_UnitTests_LL
    SOURCES
      StridedMap_UnitTests
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    DEFINES -DXPETRA_TEST_USE_LONGLONG_GO
  )

  TRIBITS_ADD_EXECUTABLE(
    StridedMapFactory_UnitTests_LL
    SOURCES
      StridedMapFactory_UnitTests
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    DEFINES -DXPETRA_TEST_USE_LONGLONG_GO
  )

  IF (${PACKAGE_NAME}_ENABLE_Epetra)

    TRIBITS_ADD_TEST(
      StridedMap_UnitTests_LL
      NAME "StridedMapUnitTestsEpetraLL"
      ARGS "--linAlgebra=Epetra"
      COMM serial mpi
      STANDARD_PASS_OUTPUT
    )

    TRIBITS_ADD_TEST(
      StridedMapFactory_UnitTests_LL
      NAME "StridedMapFactoryUnitTestsEpetraLL"
      ARGS "--linAlgebra=Epetra"
      COMM serial mpi
      STANDARD_PASS_OUTPUT
    )

  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Tpetra)

    TRIBITS_ADD_TEST(
      StridedMap_UnitTests_LL
      NAME "StridedMapUnitTestsTpetraLL"
      ARGS "--linAlgebra=Tpetra"
      COMM serial mpi
      STANDARD_PASS_OUTPUT
    )

    TRIBITS_ADD_TEST(
      StridedMapFactory_UnitTests_LL
      NAME "StridedMapFactoryUnitTestsTpetraLL"
      ARGS "--linAlgebra=Tpetra"
      COMM serial mpi
      STANDARD_PASS_OUTPUT
    )

  ENDIF()
ENDIF()

