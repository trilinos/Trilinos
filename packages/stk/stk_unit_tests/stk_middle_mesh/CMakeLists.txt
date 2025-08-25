FILE(GLOB SOURCES CONFIGURE_DEPENDS *.cpp util/*.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE(stk_middle_mesh_utest
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    COMM serial mpi
    NOEXEPREFIX
    )

  TRIBITS_ADD_TEST(stk_middle_mesh_utest
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    NOEXEPREFIX
    )
else()
  FILE(GLOB TEST_UTIL_SOURCES CONFIGURE_DEPENDS util/*.cpp)
  add_library(stk_middle_mesh_test_util ${TEST_UTIL_SOURCES})
  target_link_libraries(stk_middle_mesh_test_util stk_middle_mesh)
  target_link_libraries(stk_middle_mesh_test_util GTest::gtest)
  target_include_directories(stk_middle_mesh_test_util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories(stk_middle_mesh_test_util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/util)

  add_executable(stk_middle_mesh_utest ${SOURCES})
  target_link_libraries(stk_middle_mesh_utest stk_middle_mesh)
  target_link_libraries(stk_middle_mesh_utest stk_unit_test_utils)
  target_link_libraries(stk_middle_mesh_utest stk_util_parallel)
  target_link_libraries(stk_middle_mesh_utest stk_unit_main)

  FOREACH(NUMPROCS 1 2 3 4 6 8)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS})
    add_test(NAME "stk_middle_mesh_utest_np${NUMPROCS}" 
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS}  
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROCS} ${CMAKE_CURRENT_BINARY_DIR}/stk_middle_mesh_utest)
    set_property(TEST "stk_middle_mesh_utest_np${NUMPROCS}" PROPERTY PROCESSORS ${NUMPROCS})
  ENDFOREACH()

endif()

if(TARGET stk_middle_mesh_utest)
  target_include_directories(stk_middle_mesh_utest PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
  )

  INSTALL(TARGETS stk_middle_mesh_utest DESTINATION ${STK_INSTALL_BINDIR})
endif()
