kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR})
kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../test_common)

#At some point, we may wish to make this into real "performance
#tests, " in the sense that they can pass or fail.  At that point, use
#"CATEGORIES PERFORMANCE" to mark them as such. For now, we just
#build the executables for manual use, but don't run the tests.  They
#build correctly with or without MPI, but only run them with a single
#MPI process.

# add_component_subdirectory(batched)
add_component_subdirectory(blas)
# add_subdirectory(common)
# add_component_subdirectory(graph)
add_component_subdirectory(lapack)
add_component_subdirectory(ode)
add_component_subdirectory(sparse)
