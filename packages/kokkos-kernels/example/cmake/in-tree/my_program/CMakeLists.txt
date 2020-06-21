
add_executable(myTest test.cc)
target_link_libraries(myTest PUBLIC kokkos)
target_link_libraries(myTest PUBLIC kokkoskernels)


