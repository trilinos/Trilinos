FILE(GLOB SOURCES *.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE(stk_middle_mesh_doc_tests
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    COMM serial mpi
    NOEXEPREFIX
    )

  TRIBITS_ADD_TEST(stk_middle_mesh_doc_tests
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    NOEXEPREFIX
    )
else()
  add_executable(stk_middle_mesh_doc_tests ${SOURCES})
  target_link_libraries(stk_middle_mesh_doc_tests stk_middle_mesh)
  target_link_libraries(stk_middle_mesh_doc_tests stk_unit_test_utils)
  target_link_libraries(stk_middle_mesh_doc_tests stk_util_parallel)
  target_link_libraries(stk_middle_mesh_doc_tests stk_unit_main)

  add_test(NAME "stk_middle_mesh_doc_tests" COMMAND stk_middle_mesh_doc_tests)
endif()

if(TARGET stk_middle_mesh_doc_tests)
  target_include_directories(stk_middle_mesh_doc_tests PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
  )
endif()
