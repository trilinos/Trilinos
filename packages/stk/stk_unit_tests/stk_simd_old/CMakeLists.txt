
FILE(GLOB SOURCES CONFIGURE_DEPENDS *.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  TRIBITS_INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/stk_simd)
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/stk_simd)

  TRIBITS_ADD_EXECUTABLE(stk_simd_old_unit_tests
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    COMM serial mpi
    NOEXEPREFIX
    )

  TRIBITS_ADD_TEST(stk_simd_old_unit_tests
   ARGS "--gtest_filter=StkSimd.SimdAddSubtractMultDivide"
   COMM serial mpi
   PASS_REGULAR_EXPRESSION "PASS"
   FAIL_REGULAR_EXPRESSION "FAIL"
   NUM_MPI_PROCS 1
   NOEXEPREFIX
   )
else()
  add_executable(stk_simd_old_unit ${SOURCES})

  target_include_directories(stk_simd_old_unit PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
  )

  FOREACH(NUMPROCS 1)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS})
    add_test(NAME "stk_simd_old_unit_np${NUMPROCS}" 
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROCS} ${CMAKE_CURRENT_BINARY_DIR}/stk_simd_old_unit)
    set_property(TEST "stk_simd_old_unit_np${NUMPROCS}" PROPERTY PROCESSORS ${NUMPROCS})
  ENDFOREACH()

  target_link_libraries(stk_simd_old_unit stk_simd)
  target_link_libraries(stk_simd_old_unit stk_math)
  target_link_libraries(stk_simd_old_unit stk_unit_test_utils)
  target_link_libraries(stk_simd_old_unit stk_unit_main)

  INSTALL(TARGETS stk_simd_old_unit DESTINATION ${STK_INSTALL_BINDIR})

endif()

