# Adding source directory to the build
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/common/src)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/common/impl)

LIST(APPEND SOURCES common/src/KokkosKernels_EagerInitialize.cpp)
