kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../test_common)

kokkoskernels_add_executable_and_test(wiki_blas1_dot         SOURCES KokkosBlas1_wiki_dot.cpp)
kokkoskernels_add_executable_and_test(wiki_blas1_iamax       SOURCES KokkosBlas1_wiki_iamax.cpp)
kokkoskernels_add_executable_and_test(wiki_blas1_nrm1        SOURCES KokkosBlas1_wiki_nrm1.cpp)
kokkoskernels_add_executable_and_test(wiki_blas1_nrm2        SOURCES KokkosBlas1_wiki_nrm2.cpp)
kokkoskernels_add_executable_and_test(wiki_blas1_rotg_rot    SOURCES KokkosBlas1_wiki_rotg_rot.cpp)
kokkoskernels_add_executable_and_test(wiki_blas1_rotmg_rotm  SOURCES KokkosBlas1_wiki_rotmg_rotm.cpp)
kokkoskernels_add_executable_and_test(wiki_blas2_ger         SOURCES KokkosBlas2_wiki_ger.cpp)
kokkoskernels_add_executable_and_test(wiki_blas2_syr         SOURCES KokkosBlas2_wiki_syr.cpp)
kokkoskernels_add_executable_and_test(wiki_blas2_syr2        SOURCES KokkosBlas2_wiki_syr2.cpp)
