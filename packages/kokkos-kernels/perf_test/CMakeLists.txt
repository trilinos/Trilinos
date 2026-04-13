if(KokkosKernels_ENABLE_PERFTESTS)
  kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
  kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

  kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../test_common)

  # At some point, we may wish to make this into real "performance
  # tests, " in the sense that they can pass or fail.  At that point, use
  # "CATEGORIES PERFORMANCE" to mark them as such.For now, we just
  # build the executables for manual use, but don't run the tests.

  kokkoskernels_include_directories(sparse)

  add_component_subdirectory(batched)
  add_component_subdirectory(graph)
  add_component_subdirectory(sparse)
  add_component_subdirectory(blas)
  add_subdirectory(performance)
  #add_subdirectory(common)

endif()
