TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  02_Use_ModelEvaluator
  SOURCES 02_Use_ModelEvaluator.cpp VanDerPol_ModelEvaluator_02.hpp VanDerPol_ModelEvaluator_02.cpp
  NUM_MPI_PROCS 1
  )
