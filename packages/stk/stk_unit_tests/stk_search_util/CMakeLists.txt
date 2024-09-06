
FILE(GLOB SOURCES *.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_search_util_utest
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    ARGS "--gtest_filter=-SearchUtil.Intrepid2*"
    COMM serial mpi
    NUM_MPI_PROCS 1
  )
else()
  add_executable(stk_search_util_utest ${SOURCES})
  target_link_libraries(stk_search_util_utest PUBLIC stk_search_util)
  target_link_libraries(stk_search_util_utest PUBLIC stk_unit_test_utils)
  target_link_libraries(stk_search_util_utest PUBLIC stk_mesh_fixtures)
  target_link_libraries(stk_search_util_utest PUBLIC stk_unit_main)
endif()

