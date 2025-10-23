STK_SUBPACKAGE(STKPerformance_tests)

if(HAVE_STK_Trilinos)

  message("STK performance tests are not currently supported for Trilinos.")

else()

  FILE(GLOB SOURCES_MIDDLE_MESH CONFIGURE_DEPENDS stk_middle_mesh/*.cpp)

  FILE(GLOB SOURCES CONFIGURE_DEPENDS stk_balance/*.cpp
                                      stk_io/*.cpp
                                      stk_mesh/*.cpp
                                      stk_mesh/entity_sorting/*.cpp
                                      general/*.cpp
                                      stk_search/*.cpp
                                      stk_transfer/*.cpp
                                      stk_util/*.cpp
                                      stk_tools/*.cpp
      )

  if (STK_ENABLE_STKMiddle_mesh_util)
    LIST(APPEND SOURCES ${SOURCES_MIDDLE_MESH})
  endif()

  FILE(GLOB HEADERS CONFIGURE_DEPENDS stk_mesh/*.hpp)

  add_executable(stk_performance_tests ${SOURCES})
  target_link_libraries(stk_performance_tests stk_unit_main)
  target_link_libraries(stk_performance_tests stk_balance_lib)
  if (STK_ENABLE_STKMiddle_mesh_util)
    target_link_libraries(stk_performance_tests stk_middle_mesh_util)
    target_link_libraries(stk_performance_tests stk_middle_mesh)
  endif()
  target_link_libraries(stk_performance_tests stk_coupling)
  target_link_libraries(stk_performance_tests stk_search)
  target_link_libraries(stk_performance_tests stk_search_util)  
  target_link_libraries(stk_performance_tests stk_transfer_util)
  target_link_libraries(stk_performance_tests stk_io_util)
  target_link_libraries(stk_performance_tests stk_util_diag)
  target_link_libraries(stk_performance_tests stk_tools_lib)
  target_link_libraries(stk_performance_tests stk_mesh_fixtures)
  target_link_libraries(stk_performance_tests stk_unit_test_utils)

  target_include_directories(stk_performance_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
  
  INSTALL(TARGETS stk_performance_tests DESTINATION ${STK_INSTALL_BINDIR})

endif()

STK_SUBPACKAGE_POSTPROCESS()

