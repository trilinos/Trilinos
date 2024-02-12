# CMAKE File for "MyApp" application building against an installed Trilinos

cmake_minimum_required(VERSION 3.0)

# Define the project and the compilers
#
# NOTE: You can't call find_package(Trilinos) for a CUDA build without first
# defining the compilers.
#
project(MyApp C CXX)

# Disable Kokkos warning about not supporting C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)

# Get just Tpetra as
find_package(Tpetra REQUIRED)

# Echo trilinos build info just for fun
message("\nFound Tpetra!  Here are the details: ")
message("   Tpetra_DIR = ${Tpetra_DIR}")

# Build the APP and link to Trilinos
add_executable(MyApp ${CMAKE_CURRENT_SOURCE_DIR}/app.cpp)
target_link_libraries(MyApp  Tpetra::all_libs)

# Set up a test
enable_testing()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.xml
  ${CMAKE_CURRENT_BINARY_DIR}/input.xml COPYONLY)
set(NUM_MPI_PROCS 4)
add_test(MyAppTest mpiexec -np ${NUM_MPI_PROCS} "${CMAKE_CURRENT_BINARY_DIR}/MyApp")
set_tests_properties(MyAppTest PROPERTIES
  PROCESSORS ${NUM_MPI_PROCS}
  PASS_REGULAR_EXPRESSION "vec.norm1[(][)] = 40"
  )
# NOTE: Above, mpiexec with mpich-3.2 requires you pass in the abs path to
# MyApp or mpiexec says it can't find it, even though it is running in the
# correct directory (see #10813).
