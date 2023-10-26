TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  03_Intro_SolutionState
  SOURCES 03_Intro_SolutionState.cpp ../02_Use_ModelEvaluator/VanDerPol_ModelEvaluator_02.hpp ../02_Use_ModelEvaluator/VanDerPol_ModelEvaluator_02.cpp
  NUM_MPI_PROCS 1
  )
