INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
                    ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExplicitRK
  SOURCES Tempus_ExplicitRKTest.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExplicitRK_Combined_FSA
  SOURCES Tempus_ExplicitRK_Combined_FSA.cpp Tempus_ExplicitRK_FSA.hpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExplicitRK_Staggered_FSA
  SOURCES Tempus_ExplicitRK_Staggered_FSA.cpp Tempus_ExplicitRK_FSA.hpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExplicitRK_PseudoTransient_SA
  SOURCES Tempus_ExplicitRK_PseudoTransient_SA.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExplicitRK_ASA
  SOURCES Tempus_ExplicitRK_ASA.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Test_ExplicitRK_CopyFiles
  DEST_FILES Tempus_ExplicitRK_SinCos.xml Tempus_ExplicitRK_SteadyQuadratic.xml
  EXEDEPS ExplicitRK
  )
