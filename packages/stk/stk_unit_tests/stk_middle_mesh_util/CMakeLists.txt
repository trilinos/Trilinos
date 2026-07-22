FILE(GLOB SOURCES CONFIGURE_DEPENDS *.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE(stk_middle_mesh_util_utest
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    COMM serial mpi
    NOEXEPREFIX
    )

  TRIBITS_ADD_TEST(stk_middle_mesh_util_utest
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    NOEXEPREFIX
    )
else()
  add_executable(stk_middle_mesh_util_utest ${SOURCES})
  target_link_libraries(stk_middle_mesh_util_utest stk_middle_mesh_util)
  target_link_libraries(stk_middle_mesh_util_utest stk_middle_mesh)
  target_link_libraries(stk_middle_mesh_util_utest stk_mesh_base)
  target_link_libraries(stk_middle_mesh_util_utest stk_topology)
  target_link_libraries(stk_middle_mesh_util_utest stk_unit_test_utils)
  target_link_libraries(stk_middle_mesh_util_utest stk_unit_main)

  FOREACH(NUMPROCS 1 2)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS})
    add_test(NAME "stk_middle_mesh_util_utest_np${NUMPROCS}" 
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS} 
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROCS} ${CMAKE_CURRENT_BINARY_DIR}/stk_middle_mesh_util_utest 
             --gtest_filter=-*externalFile)
    set_property(TEST "stk_middle_mesh_util_utest_np${NUMPROCS}" PROPERTY PROCESSORS ${NUMPROCS})
  ENDFOREACH()

endif()

if(TARGET stk_middle_mesh_util_utest)
  target_include_directories(stk_middle_mesh_util_utest PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
  )

  INSTALL(TARGETS stk_middle_mesh_util_utest DESTINATION ${STK_INSTALL_BINDIR})
endif()
