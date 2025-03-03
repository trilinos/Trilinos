
FILE(GLOB SOURCES_A CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/utest_a/*.cpp)
FILE(GLOB SOURCES_B CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/utest_b/*.cpp)
FILE(GLOB SOURCES_C CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/utest_c/*.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_topology_utest_a
    SOURCES ${SOURCES_A}
    TESTONLYLIBS stk_unit_main
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    NOEXEPREFIX
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_topology_utest_b
    SOURCES ${SOURCES_B}
    TESTONLYLIBS stk_unit_main
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    NOEXEPREFIX
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_topology_utest_c
    SOURCES ${SOURCES_C}
    TESTONLYLIBS stk_unit_main
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    NOEXEPREFIX
    )
else()
  add_executable(stk_topology_utest_a ${SOURCES_A})
  target_link_libraries(stk_topology_utest_a stk_topology)
  target_link_libraries(stk_topology_utest_a stk_unit_test_utils)
  target_link_libraries(stk_topology_utest_a stk_util_env)
  target_link_libraries(stk_topology_utest_a stk_util_ngp)
  target_link_libraries(stk_topology_utest_a stk_util_util)
  target_link_libraries(stk_topology_utest_a stk_unit_main)

  FOREACH(NUMPROCS 1)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/a_np${NUMPROCS})
    add_test(NAME "stk_topology_utest_a_np${NUMPROCS}" 
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/a_np${NUMPROCS}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROCS} ${CMAKE_CURRENT_BINARY_DIR}/stk_topology_utest_a)
    set_property(TEST "stk_topology_utest_a_np${NUMPROCS}" PROPERTY PROCESSORS ${NUMPROCS})
  ENDFOREACH()

  add_executable(stk_topology_utest_b ${SOURCES_B})
  target_link_libraries(stk_topology_utest_b stk_topology)
  target_link_libraries(stk_topology_utest_b stk_unit_test_utils)
  target_link_libraries(stk_topology_utest_b stk_util_env)
  target_link_libraries(stk_topology_utest_b stk_util_ngp)
  target_link_libraries(stk_topology_utest_b stk_util_util)
  target_link_libraries(stk_topology_utest_b stk_unit_main)

  FOREACH(NUMPROCS 1)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/b_np${NUMPROCS})
    add_test(NAME "stk_topology_utest_b_np${NUMPROCS}" 
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/b_np${NUMPROCS}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROCS} ${CMAKE_CURRENT_BINARY_DIR}/stk_topology_utest_b)
    set_property(TEST "stk_topology_utest_b_np${NUMPROCS}" PROPERTY PROCESSORS ${NUMPROCS})
  ENDFOREACH()

  add_executable(stk_topology_utest_c ${SOURCES_C})
  target_link_libraries(stk_topology_utest_c stk_topology)
  target_link_libraries(stk_topology_utest_c stk_unit_test_utils)
  target_link_libraries(stk_topology_utest_c stk_util_env)
  target_link_libraries(stk_topology_utest_c stk_util_ngp)
  target_link_libraries(stk_topology_utest_c stk_util_util)
  target_link_libraries(stk_topology_utest_c stk_unit_main)

  FOREACH(NUMPROCS 1)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/c_np${NUMPROCS})
    add_test(NAME "stk_topology_utest_c_np${NUMPROCS}" 
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/c_np${NUMPROCS}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROCS} ${CMAKE_CURRENT_BINARY_DIR}/stk_topology_utest_c)
    set_property(TEST "stk_topology_utest_c_np${NUMPROCS}" PROPERTY PROCESSORS ${NUMPROCS})
  ENDFOREACH()

endif()

if(TARGET stk_topology_utest_a)
target_include_directories(stk_topology_utest_a PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

  INSTALL(TARGETS stk_topology_utest_a DESTINATION ${STK_INSTALL_BINDIR})
endif()

if(TARGET stk_topology_utest_b)
target_include_directories(stk_topology_utest_b PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

  INSTALL(TARGETS stk_topology_utest_b DESTINATION ${STK_INSTALL_BINDIR})
endif()

if(TARGET stk_topology_utest_c)
target_include_directories(stk_topology_utest_c PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

  INSTALL(TARGETS stk_topology_utest_c DESTINATION ${STK_INSTALL_BINDIR})
endif()
