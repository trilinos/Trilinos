
FILE(GLOB SOURCES *.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_transfer_util_unit_tests
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
  )
else()
  add_executable(stk_transfer_util_unit_tests ${SOURCES})
  target_link_libraries(stk_transfer_util_unit_tests PUBLIC stk_transfer_util)
  target_link_libraries(stk_transfer_util_unit_tests PUBLIC stk_unit_test_utils)
  target_link_libraries(stk_transfer_util_unit_tests PUBLIC stk_io)
  target_link_libraries(stk_transfer_util_unit_tests PUBLIC stk_unit_main)
endif()

