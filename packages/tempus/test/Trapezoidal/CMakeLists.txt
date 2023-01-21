TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
                    ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Trapezoidal
  SOURCES Tempus_TrapezoidalTest.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Test_Trapezoidal_CopyFiles
  DEST_FILES Tempus_Trapezoidal_SinCos.xml Tempus_Trapezoidal_VanDerPol.xml
  EXEDEPS Trapezoidal
  )
