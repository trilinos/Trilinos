kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../test_common)

kokkoskernels_add_executable_and_test(wiki_coloring   SOURCES KokkosGraph_wiki_coloring.cpp)
kokkoskernels_add_executable_and_test(wiki_mis2       SOURCES KokkosGraph_wiki_mis2.cpp)
kokkoskernels_add_executable_and_test(wiki_coarsening SOURCES KokkosGraph_wiki_coarsening.cpp)
kokkoskernels_add_executable_and_test(wiki_rcm        SOURCES KokkosGraph_wiki_rcm.cpp)
