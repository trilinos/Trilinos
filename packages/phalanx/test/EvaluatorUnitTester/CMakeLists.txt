TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  evaluator_unit_tester
  SOURCES ${PHALANX_UNIT_TEST_MAIN}
          EvaluatorUnitTester.cpp
          SimpleEvaluator.hpp
          SimpleEvaluator_Def.hpp
          DuplicateFieldEvaluator.hpp
          DuplicateFieldEvaluator_Def.hpp
          MyTraits.hpp
  COMM serial mpi
  NUM_MPI_PROCS 1
)

