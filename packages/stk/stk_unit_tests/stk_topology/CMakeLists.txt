
FILE(GLOB SOURCES *.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_topology_unit_tests
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    NOEXEPREFIX
    )
else()
  add_executable(stk_topology_unit_tests ${SOURCES})
  target_link_libraries(stk_topology_unit_tests stk_topology)
  target_link_libraries(stk_topology_unit_tests stk_unit_test_utils)
  target_link_libraries(stk_topology_unit_tests stk_util_env)
  target_link_libraries(stk_topology_unit_tests stk_util_ngp)
  target_link_libraries(stk_topology_unit_tests stk_util_util)
  target_link_libraries(stk_topology_unit_tests stk_unit_main)

  add_test(NAME "stk_topology_unit_tests" COMMAND stk_topology_unit_tests)
endif()

target_include_directories(stk_topology_unit_tests PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

