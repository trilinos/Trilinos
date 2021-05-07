# CMAKE File for "MyApp" application building against an installed Trilinos

cmake_minimum_required(VERSION 3.17.1)

# Declare project but don't process compilers yet
#
# You must do this before calling find_package(Trilinos) when
# BUILD_SHARED_LIBS=ON or you get a lot of nasty warnings.
#
project(MyApp NONE)

# Disable Kokkos warning about not supporting C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)

# Get Trilinos as one entity but require the packages being used
find_package(Trilinos REQUIRED COMPONENTS Teuchos Tpetra)

# Echo trilinos build info just for fun
MESSAGE("\nFound Trilinos!  Here are the details: ")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("End of Trilinos details\n")

# Make sure to use same compilers and flags as Trilinos
set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

set(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
set(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
set(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")

# End of setup and error checking.

# Now enable the compilers now that we have gotten them from Trilinos
enable_language(C)
enable_language(CXX)
if (CMAKE_Fortran_COMPILER)
  enable_language(Fortran)
endif()
# NOTE: Technically this project only has C++ code so we could only call
# enable_language(CXX) here.  But we leave the enable of C and Fortran here to
# demonstrate how to do this.  Also, note that Fortran is optional in Trilinos
# so we make the enable of Fortran optional based on it's enable in Trilinos.

# Build the APP and link to Trilinos
add_executable(MyApp ${CMAKE_CURRENT_SOURCE_DIR}/app.cpp)
target_include_directories(MyApp PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR} ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
target_link_libraries(MyApp ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}) 

# Set up a test
enable_testing()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.xml
  ${CMAKE_CURRENT_BINARY_DIR}/input.xml COPYONLY)
set(NUM_MPI_PROCS 4)
add_test(MyAppTest mpiexec -np ${NUM_MPI_PROCS} MyApp)
set_tests_properties(MyAppTest PROPERTIES
  PROCESSORS ${NUM_MPI_PROCS}
  PASS_REGULAR_EXPRESSION "vec.norm1[(][)] = 40"
  )
