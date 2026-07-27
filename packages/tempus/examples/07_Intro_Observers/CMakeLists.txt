TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  07_Intro_Observers
  SOURCES 07_Intro_Observers.cpp TutorialObserver.hpp ../02_Use_ModelEvaluator/VanDerPol_ModelEvaluator_02.hpp ../02_Use_ModelEvaluator/VanDerPol_ModelEvaluator_02.cpp
  NUM_MPI_PROCS 1
  )
