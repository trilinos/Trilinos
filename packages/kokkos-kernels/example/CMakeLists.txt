kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../test_common)

#ADD_SUBDIRECTORY(graph)
add_subdirectory(wiki)
add_subdirectory(gmres)
add_subdirectory(batched_solve)
add_subdirectory(half)
