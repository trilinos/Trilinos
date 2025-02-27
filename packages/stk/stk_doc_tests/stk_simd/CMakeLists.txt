
FILE(GLOB SOURCES CONFIGURE_DEPENDS *.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE(stk_simd_docs_tests
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    COMM serial mpi
    NOEXEPREFIX
  )

  TRIBITS_ADD_TEST(
    stk_simd_docs_tests
    ARGS ""
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
    NUM_MPI_PROCS 1
    NOEXEPREFIX
  )
else()
  add_executable(stk_simd_docs_tests ${SOURCES})

  target_link_libraries(stk_simd_docs_tests stk_simd)
  target_link_libraries(stk_simd_docs_tests stk_unit_test_utils)
  target_link_libraries(stk_simd_docs_tests stk_unit_main)
  
  FOREACH(NUMPROCS 1)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS})
    add_test(NAME "stk_simd_docs_tests_np${NUMPROCS}"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS} 
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROCS} ${CMAKE_CURRENT_BINARY_DIR}/stk_simd_docs_tests)
    set_property(TEST "stk_simd_docs_tests_np${NUMPROCS}" PROPERTY PROCESSORS ${NUMPROCS})
  ENDFOREACH()

  INSTALL(TARGETS stk_simd_docs_tests DESTINATION ${STK_INSTALL_BINDIR})

endif()

if(TARGET stk_simd_docs_tests)
  target_include_directories(stk_simd_docs_tests PRIVATE
    $<BUILD_INTERFACE:${STK_TOPLEVEL_SOURCE_DIR}/stk_simd>
    $<INSTALL_INTERFACE:include>
  )

endif()
