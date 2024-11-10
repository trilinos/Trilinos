
FILE(GLOB SOURCES *.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE(stk_simd_doc_test
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    COMM serial mpi
    NOEXEPREFIX
  )

  TRIBITS_ADD_TEST(
    stk_simd_doc_test
    ARGS ""
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
    NUM_MPI_PROCS 1
    NOEXEPREFIX
  )
else()
  add_executable(stk_simd_doc_test ${SOURCES})

  target_link_libraries(stk_simd_doc_test stk_simd)
  target_link_libraries(stk_simd_doc_test stk_unit_test_utils)
  target_link_libraries(stk_simd_doc_test stk_unit_main)
endif()

if(TARGET stk_simd_doc_test)
  target_include_directories(stk_simd_doc_test PRIVATE
    $<BUILD_INTERFACE:${STK_TOPLEVEL_SOURCE_DIR}/stk_simd>
    $<INSTALL_INTERFACE:include>
  )
endif()
