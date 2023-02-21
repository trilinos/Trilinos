TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
                    ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ForwardEuler
  SOURCES Tempus_ForwardEulerTest.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Test_ForwardEuler_CopyFiles
  DEST_FILES Tempus_ForwardEuler_SinCos.xml Tempus_ForwardEuler_VanDerPol.xml Tempus_ForwardEuler_NumberOfTimeSteps.xml
  EXEDEPS ForwardEuler
  )
