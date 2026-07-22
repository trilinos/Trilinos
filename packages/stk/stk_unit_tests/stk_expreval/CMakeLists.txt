
FILE(GLOB SOURCES CONFIGURE_DEPENDS *.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_expreval_unit_tests
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
  )
else()
  add_executable(stk_expreval_utest ${SOURCES})
  target_link_libraries(stk_expreval_utest stk_expreval)
  target_link_libraries(stk_expreval_utest stk_ngp_test)
  target_link_libraries(stk_expreval_utest stk_unit_main)

  INSTALL(TARGETS stk_expreval_utest DESTINATION ${STK_INSTALL_BINDIR})

  FOREACH(NUMPROCS 1)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS})
    add_test(NAME "stk_expreval_utest_np${NUMPROCS}"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROCS} ${CMAKE_CURRENT_BINARY_DIR}/stk_expreval_utest)
    set_property(TEST "stk_expreval_utest_np${NUMPROCS}" PROPERTY PROCESSORS ${NUMPROCS})
  ENDFOREACH()
  
endif()

