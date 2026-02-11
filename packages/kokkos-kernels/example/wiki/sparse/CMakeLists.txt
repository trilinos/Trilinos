kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../test_common)

if(KOKKOSKERNELS_ENABLE_EXPERIMENTAL)
  kokkoskernels_add_executable_and_test(wiki_bsrmatrix SOURCES KokkosSparse_wiki_bsrmatrix.cpp)
endif()

kokkoskernels_add_executable_and_test(wiki_bsrmatrix_2   SOURCES KokkosSparse_wiki_bsrmatrix_2.cpp)
kokkoskernels_add_executable_and_test(wiki_crsmatrix     SOURCES KokkosSparse_wiki_crsmatrix.cpp)
kokkoskernels_add_executable_and_test(wiki_spmv          SOURCES KokkosSparse_wiki_spmv.cpp)
kokkoskernels_add_executable_and_test(wiki_spadd         SOURCES KokkosSparse_wiki_spadd.cpp)
kokkoskernels_add_executable_and_test(wiki_spgemm        SOURCES KokkosSparse_wiki_spgemm.cpp)
kokkoskernels_add_executable_and_test(wiki_gauss_seidel  SOURCES KokkosSparse_wiki_gauss_seidel.cpp)
