# Adding source directory to the build
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/common/src)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/common/impl)

list(APPEND SOURCES common/src/KokkosKernels_EagerInitialize.cpp)
