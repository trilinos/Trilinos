
FILE(GLOB SOURCES_A ${CMAKE_CURRENT_SOURCE_DIR}/utest_a/*.cpp)
FILE(GLOB SOURCES_B ${CMAKE_CURRENT_SOURCE_DIR}/utest_b/*.cpp)
FILE(GLOB SOURCES_C ${CMAKE_CURRENT_SOURCE_DIR}/utest_c/*.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_topology_unit_tests_a
    SOURCES ${SOURCES_A}
    TESTONLYLIBS stk_unit_main
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    NOEXEPREFIX
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_topology_unit_tests_b
    SOURCES ${SOURCES_B}
    TESTONLYLIBS stk_unit_main
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    NOEXEPREFIX
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_topology_unit_tests_c
    SOURCES ${SOURCES_C}
    TESTONLYLIBS stk_unit_main
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    NOEXEPREFIX
    )
else()
  add_executable(stk_topology_unit_tests_a ${SOURCES_A})
  target_link_libraries(stk_topology_unit_tests_a stk_topology)
  target_link_libraries(stk_topology_unit_tests_a stk_unit_test_utils)
  target_link_libraries(stk_topology_unit_tests_a stk_util_env)
  target_link_libraries(stk_topology_unit_tests_a stk_util_ngp)
  target_link_libraries(stk_topology_unit_tests_a stk_util_util)
  target_link_libraries(stk_topology_unit_tests_a stk_unit_main)

  add_test(NAME "stk_topology_unit_tests_a" COMMAND stk_topology_unit_tests_a)

  add_executable(stk_topology_unit_tests_b ${SOURCES_B})
  target_link_libraries(stk_topology_unit_tests_b stk_topology)
  target_link_libraries(stk_topology_unit_tests_b stk_unit_test_utils)
  target_link_libraries(stk_topology_unit_tests_b stk_util_env)
  target_link_libraries(stk_topology_unit_tests_b stk_util_ngp)
  target_link_libraries(stk_topology_unit_tests_b stk_util_util)
  target_link_libraries(stk_topology_unit_tests_b stk_unit_main)

  add_test(NAME "stk_topology_unit_tests_b" COMMAND stk_topology_unit_tests_b)

  add_executable(stk_topology_unit_tests_c ${SOURCES_C})
  target_link_libraries(stk_topology_unit_tests_c stk_topology)
  target_link_libraries(stk_topology_unit_tests_c stk_unit_test_utils)
  target_link_libraries(stk_topology_unit_tests_c stk_util_env)
  target_link_libraries(stk_topology_unit_tests_c stk_util_ngp)
  target_link_libraries(stk_topology_unit_tests_c stk_util_util)
  target_link_libraries(stk_topology_unit_tests_c stk_unit_main)

  add_test(NAME "stk_topology_unit_tests_c" COMMAND stk_topology_unit_tests_c)
endif()

if(TARGET stk_topology_unit_tests_a)
target_include_directories(stk_topology_unit_tests_a PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
endif()

if(TARGET stk_topology_unit_tests_b)
target_include_directories(stk_topology_unit_tests_b PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
endif()

if(TARGET stk_topology_unit_tests_c)
target_include_directories(stk_topology_unit_tests_c PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
endif()
