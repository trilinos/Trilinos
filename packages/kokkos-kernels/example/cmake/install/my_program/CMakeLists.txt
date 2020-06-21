
add_executable(myTest test.cc)

# I now link to the CMake-specific target
target_link_libraries(myTest PUBLIC Kokkos::kokkoskernels)


