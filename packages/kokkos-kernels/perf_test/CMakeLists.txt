INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../test_common)

# At some point, we may wish to make this into real "performance
# tests," in the sense that they can pass or fail.  At that point, use
# "CATEGORIES PERFORMANCE" to mark them as such.  For now, we just
# build the executables for manual use, but don't run the tests.  They
# build correctly with or without MPI, but only run them with a single
# MPI process.

ADD_SUBDIRECTORIES(graph)
ADD_SUBDIRECTORIES(sparse)
#ADD_SUBDIRECTORIES(common)

