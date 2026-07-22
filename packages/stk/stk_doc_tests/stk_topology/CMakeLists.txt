
FILE(GLOB SOURCES CONFIGURE_DEPENDS *.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_topology_doc_test
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    COMM serial mpi
    NUM_MPI_PROCS 1
    NOEXEPREFIX
  )
else()
  add_executable(stk_topology_doc_test ${SOURCES})
  target_link_libraries(stk_topology_doc_test stk_topology)
  target_link_libraries(stk_topology_doc_test stk_unit_test_utils)
  target_link_libraries(stk_topology_doc_test stk_util_env)
  target_link_libraries(stk_topology_doc_test stk_util_ngp)
  target_link_libraries(stk_topology_doc_test stk_util_util)
  target_link_libraries(stk_topology_doc_test stk_unit_main)

  FOREACH(NUMPROCS 1)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS})
    add_test(NAME "stk_topology_doc_test_np${NUMPROCS}"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS} 
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROCS} ${CMAKE_CURRENT_BINARY_DIR}/stk_topology_doc_test)
    set_property(TEST "stk_topology_doc_test_np${NUMPROCS}" PROPERTY PROCESSORS ${NUMPROCS})
  ENDFOREACH()

  INSTALL(TARGETS stk_topology_doc_test DESTINATION ${STK_INSTALL_BINDIR})
endif()

