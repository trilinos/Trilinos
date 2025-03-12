FILE(GLOB SOURCES CONFIGURE_DEPENDS *.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_emend_unit_tests
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    NOEXEPREFIX
    )
else()
  add_executable(stk_emend_independent_set_utest ${SOURCES})
  target_link_libraries(stk_emend_independent_set_utest stk_emend)
  target_link_libraries(stk_emend_independent_set_utest stk_util_util)
  target_link_libraries(stk_emend_independent_set_utest stk_util_parallel)
  target_link_libraries(stk_emend_independent_set_utest stk_unit_test_utils)
  target_link_libraries(stk_emend_independent_set_utest stk_unit_main)

  target_include_directories(stk_emend_independent_set_utest PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
  )

  INSTALL(TARGETS stk_emend_independent_set_utest DESTINATION ${STK_INSTALL_BINDIR})

  FOREACH(NUMPROCS 1 2 3 4)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS})
    add_test(NAME "stk_emend_independent_set_utest_np${NUMPROCS}" 
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROCS} ${CMAKE_CURRENT_BINARY_DIR}/stk_emend_independent_set_utest)
    set_property(TEST "stk_emend_independent_set_utest_np${NUMPROCS}" PROPERTY PROCESSORS ${NUMPROCS})
  ENDFOREACH()
  
endif()
