INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
                    ${CMAKE_CURRENT_SOURCE_DIR})

#IKT, replace TRILINOS_ADD_EXECUTABLE code with 
#the following, once test is running properly. 

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Newmark
  SOURCES Tempus_NewmarkTest.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

#TRIBITS_ADD_EXECUTABLE(
#  Newmark
#  SOURCES Tempus_NewmarkTest.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN} 
#  TESTONLYLIBS tempus_test_models
#  COMM serial mpi 
#  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Test_NewmarkImplicitAForm_CopyFiles
  DEST_FILES Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_SecondOrder.xml Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_FirstOrder.xml Tempus_NewmarkExplicitAForm_SinCos.xml Tempus_NewmarkExplicitAForm_BallParabolic.xml Tempus_NewmarkExplicitAForm_HarmonicOscillator_Damped.xml Tempus_NewmarkImplicitDForm_HarmonicOscillator_Damped_SecondOrder.xml 
  EXEDEPS Newmark
  )
