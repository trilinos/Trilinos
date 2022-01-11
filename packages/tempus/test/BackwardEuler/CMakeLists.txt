INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
                    ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BackwardEuler
  SOURCES Tempus_BackwardEulerTest.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BackwardEuler_Combined_FSA
  SOURCES Tempus_BackwardEuler_Combined_FSA.cpp Tempus_BackwardEuler_FSA.hpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BackwardEuler_Staggered_FSA
  SOURCES Tempus_BackwardEuler_Staggered_FSA.cpp Tempus_BackwardEuler_FSA.hpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BackwardEuler_PseudoTransient_SA
  SOURCES Tempus_BackwardEuler_PseudoTransient_SA.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BackwardEuler_ASA
  SOURCES Tempus_BackwardEuler_ASA.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Test_BackwardEuler_CopyFiles
  DEST_FILES Tempus_BackwardEuler_SinCos.xml Tempus_BackwardEuler_SinCos_ASA.xml Tempus_BackwardEuler_CDR.xml Tempus_BackwardEuler_VanDerPol.xml Tempus_BackwardEuler_SteadyQuadratic.xml
  EXEDEPS BackwardEuler
  )
