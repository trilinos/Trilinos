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

# Get Trilinos as one entity but require just the packages we are using
find_package(Trilinos REQUIRED COMPONENTS Teuchos Tpetra)

# Echo trilinos build info just for fun
MESSAGE("\nFound Trilinos!  Here are the details: ")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("End of Trilinos details\n")

# Build the APP and link to Trilinos
add_executable(MyApp ${CMAKE_CURRENT_SOURCE_DIR}/app.cpp)
target_link_libraries(MyApp  Trilinos::all_selected_libs)
# Or, above could have linked to just Tpetra::all_libs

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
