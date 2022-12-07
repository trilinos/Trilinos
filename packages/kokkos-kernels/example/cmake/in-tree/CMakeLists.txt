#Kokkos requires at least 3.16
cmake_minimum_required (VERSION 3.16)

project (MyProgram)

# The add_subdirectory calls must occur in the order below

# This assumes you have copied the source, sym-linked to source,
# or have kokkos available as a git submodule
add_subdirectory(kokkos)

# This assumes you have copied the source, sym-linked to source,
# or have kokkos-kernels available as a git submodule
add_subdirectory(kokkos-kernels)

# Now I can add code for my specific program
add_subdirectory(my_program)

