
FILE(GLOB SOURCES *.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_expreval_unit_tests
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
  )
else()
  add_executable(stk_expreval_unit_tests ${SOURCES})
  target_link_libraries(stk_expreval_unit_tests stk_expreval)
  target_link_libraries(stk_expreval_unit_tests stk_ngp_test)
  target_link_libraries(stk_expreval_unit_tests stk_unit_main)

  add_test(NAME "stk_expreval_unit_tests" COMMAND stk_expreval_unit_tests)
endif()

