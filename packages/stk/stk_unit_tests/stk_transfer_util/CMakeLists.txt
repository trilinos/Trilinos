
FILE(GLOB SOURCES CONFIGURE_DEPENDS *.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_transfer_util_unit_tests
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
  )
else()
  add_executable(stk_transfer_util_utest ${SOURCES})
  target_link_libraries(stk_transfer_util_utest PUBLIC stk_transfer_util)
  target_link_libraries(stk_transfer_util_utest PUBLIC stk_unit_test_utils)
  target_link_libraries(stk_transfer_util_utest PUBLIC stk_io)
  target_link_libraries(stk_transfer_util_utest PUBLIC stk_unit_main)

  FOREACH(NUMPROCS 1 2 3 4)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS})
    add_test(NAME "stk_transfer_util_utest_np${NUMPROCS}" 
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROCS} ${CMAKE_CURRENT_BINARY_DIR}/stk_transfer_util_utest)
    set_property(TEST "stk_transfer_util_utest_np${NUMPROCS}" PROPERTY PROCESSORS ${NUMPROCS})
  ENDFOREACH()

  INSTALL(TARGETS stk_transfer_util_utest DESTINATION ${STK_INSTALL_BINDIR})
endif()

