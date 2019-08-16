INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
                    ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ForwardEulerUnitTest
  SOURCES Tempus_ForwardEulerUnitTest.cpp Tempus_UtilsUnitTest.hpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BackwardEulerUnitTest
  SOURCES Tempus_BackwardEulerUnitTest.cpp Tempus_UtilsUnitTest.hpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TrapezoidalUnitTest
  SOURCES Tempus_TrapezoidalUnitTest.cpp Tempus_UtilsUnitTest.hpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BDF2UnitTest
  SOURCES Tempus_BDF2UnitTest.cpp Tempus_UtilsUnitTest.hpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  NewmarkImplicitAFormUnitTest
  SOURCES Tempus_NewmarkImplicitAFormUnitTest.cpp Tempus_UtilsUnitTest.hpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  NewmarkImplicitDFormUnitTest
  SOURCES Tempus_NewmarkImplicitDFormUnitTest.cpp Tempus_UtilsUnitTest.hpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  NewmarkExplicitAFormUnitTest
  SOURCES Tempus_NewmarkExplicitAFormUnitTest.cpp Tempus_UtilsUnitTest.hpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  HHTAlphaUnitTest
  SOURCES Tempus_HHTAlphaUnitTest.cpp Tempus_UtilsUnitTest.hpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LeapfrogUnitTest
  SOURCES Tempus_LeapfrogUnitTest.cpp Tempus_UtilsUnitTest.hpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExplicitRKUnitTest
  SOURCES Tempus_ExplicitRKUnitTest.cpp Tempus_UtilsUnitTest.hpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DIRKUnitTest
  SOURCES Tempus_DIRKUnitTest.cpp Tempus_UtilsUnitTest.hpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  IMEX_RKUnitTest
  SOURCES Tempus_IMEX_RKUnitTest.cpp Tempus_UtilsUnitTest.hpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  IMEX_RK_PartitionedUnitTest
  SOURCES Tempus_IMEX_RK_PartitionedUnitTest.cpp Tempus_UtilsUnitTest.hpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

