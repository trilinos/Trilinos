
FILE(GLOB SOURCES CONFIGURE_DEPENDS *.cpp)

IF(NOT STK_HAS_MPI)
   LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/CommSparseHowTo.cpp")
   LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/filenameSubstitution.cpp")
   LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/usingScheduler.cpp")
ENDIF()

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_util_doc_test
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    COMM serial mpi
    NUM_MPI_PROCS 1
    NOEXEPREFIX
  )
else()
  add_executable(stk_util_doc_test ${SOURCES})
  target_link_libraries(stk_util_doc_test stk_unit_test_utils)
  target_link_libraries(stk_util_doc_test stk_util_registry)
  target_link_libraries(stk_util_doc_test stk_util_diag)
  target_link_libraries(stk_util_doc_test stk_util_env)
  target_link_libraries(stk_util_doc_test stk_util_parallel)
  target_link_libraries(stk_util_doc_test stk_util_command_line)
  target_link_libraries(stk_util_doc_test stk_unit_main)
  
  INSTALL(TARGETS stk_util_doc_test DESTINATION ${STK_INSTALL_BINDIR})

  FOREACH(NUMPROCS 1 2)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS})
    add_test(NAME "stk_util_doc_test_np${NUMPROCS}"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS} 
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROCS} ${CMAKE_CURRENT_BINARY_DIR}/stk_util_doc_test)
    set_property(TEST "stk_util_doc_test_np${NUMPROCS}" PROPERTY PROCESSORS ${NUMPROCS})
  ENDFOREACH()

endif()

