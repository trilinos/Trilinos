TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  template_manager
  SOURCES TemplateManagerTest.cpp
          MyTraits.hpp
  TESTONLYLIBS phalanx_unit_test_main
  COMM serial mpi
  NUM_MPI_PROCS 1
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  template_manager_common_api
  SOURCES TemplateManagerCommonApiTest.cpp
  TESTONLYLIBS phalanx_unit_test_main
  COMM serial mpi
  NUM_MPI_PROCS 1
)
