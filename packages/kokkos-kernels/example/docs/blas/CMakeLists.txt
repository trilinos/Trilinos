kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../test_common)

kokkoskernels_add_executable_and_test(docs_blas1_abs         SOURCES KokkosBlas1_docs_abs.cpp)
kokkoskernels_add_executable_and_test(docs_blas1_scal        SOURCES KokkosBlas1_docs_scal.cpp)
kokkoskernels_add_executable_and_test(docs_blas1_swap        SOURCES KokkosBlas1_docs_swap.cpp)
kokkoskernels_add_executable_and_test(docs_blas1_nrminf      SOURCES KokkosBlas1_docs_nrminf.cpp)
