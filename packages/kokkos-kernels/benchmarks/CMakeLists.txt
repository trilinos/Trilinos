KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../test_common)

#At some point, we may wish to make this into real "performance
#tests, " in the sense that they can pass or fail.  At that point, use
#"CATEGORIES PERFORMANCE" to mark them as such. For now, we just
#build the executables for manual use, but don't run the tests.  They
#build correctly with or without MPI, but only run them with a single
#MPI process.

# ADD_COMPONENT_SUBDIRECTORY(batched)
ADD_COMPONENT_SUBDIRECTORY(blas)
# ADD_SUBDIRECTORY(common)
# ADD_COMPONENT_SUBDIRECTORY(graph)
ADD_COMPONENT_SUBDIRECTORY(lapack)
ADD_COMPONENT_SUBDIRECTORY(ode)
ADD_COMPONENT_SUBDIRECTORY(sparse)
