#Kokkos requires at least 3.10 
#but really you should use 3.12
cmake_minimum_required (VERSION 3.10)

IF(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")
  MESSAGE(STATUS "Setting policy CMP0074 to use <Package>_ROOT variables")
  CMAKE_POLICY(SET CMP0074 NEW)
ENDIF()

project (MyProgram)

# Try to find KokkosKernels (this will also find Kokkos)
# Make sure when you run CMake you give it
# -DKokkosKernels_ROOT=<KokkosKernels_Prefix>
find_package(KokkosKernels REQUIRED)

# Now I can add code for my specific program
add_subdirectory(my_program)

