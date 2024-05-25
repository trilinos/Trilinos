
FILE(GLOB SOURCES *.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  TRIBITS_INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/stk_simd)
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/stk_simd)

  TRIBITS_ADD_EXECUTABLE(stk_simd_unit_tests
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    COMM serial mpi
    NOEXEPREFIX
    )

  TRIBITS_ADD_TEST(stk_simd_unit_tests
   COMM serial mpi
   PASS_REGULAR_EXPRESSION "PASS"
   FAIL_REGULAR_EXPRESSION "FAIL"
   NUM_MPI_PROCS 1
   NOEXEPREFIX
   )
else()
  add_executable(stk_simd_unit_tests ${SOURCES})

  target_include_directories(stk_simd_unit_tests PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
  )

  target_link_libraries(stk_simd_unit_tests stk_simd)
  target_link_libraries(stk_simd_unit_tests stk_unit_test_utils)
  target_link_libraries(stk_simd_unit_tests stk_unit_main)
endif()

