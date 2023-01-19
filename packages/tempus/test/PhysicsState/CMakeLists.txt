TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
                    ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PhysicsState
  SOURCES Tempus_PhysicsStateTest.cpp Tempus_PhysicsStateTest_StepperForwardEuler.hpp Tempus_PhysicsStateCounter.hpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Test_PhysicsState_CopyFiles
  DEST_FILES Tempus_PhysicsState_SinCos.xml
  EXEDEPS PhysicsState
  )
