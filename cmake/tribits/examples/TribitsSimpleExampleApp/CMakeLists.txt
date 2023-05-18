cmake_minimum_required(VERSION 3.23.1)

project(TribitsSimpleExApp
  DESCRIPTION
    "Example raw CMake project using packages installed from TribitsExampleProject"
  VERSION 0.0.0
  LANGUAGES NONE  # Defined below after reading in compilers
  )

find_package(TribitsExProj REQUIRED
  COMPONENTS SimpleCxx MixedLang WithSubpackages)

message("Setting compilers and flags read in from 'TribitsExProjConfig.cmake' file:")

set(CMAKE_CXX_COMPILER ${TribitsExProj_CXX_COMPILER} )
set(CMAKE_C_COMPILER ${TribitsExProj_C_COMPILER} )
set(CMAKE_Fortran_COMPILER ${TribitsExProj_Fortran_COMPILER} )

set(CMAKE_CXX_FLAGS "${TribitsExProj_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
set(CMAKE_C_FLAGS "${TribitsExProj_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
set(CMAKE_Fortran_FLAGS "${TribitsExProj_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")

# Enable the compilers now that we have gotten them from TribitsExProjConfig.cmake
enable_language(C)
enable_language(CXX)
if (CMAKE_Fortran_COMPILER)
  enable_language(Fortran)
endif()

# Build a utility and link to libraries from WithSubpackagesB
add_executable(util util.cpp)
target_link_libraries(util PRIVATE WithSubpackagesB::all_libs)

# Build the APP and link to libraries from TribitsExProj packages
add_executable(app app.cpp)
target_link_libraries(app PRIVATE TribitsExProj::all_selected_libs)

# Set up tests

enable_testing()

if (TARGET SimpleTpl::all_libs)
  set(simpleCxxDeps "simpletpl headeronlytpl")
else()
  set(simpleCxxDeps "headeronlytpl")
endif()

add_test(NAME util_test COMMAND util)
set_tests_properties(util_test PROPERTIES
  PASS_REGULAR_EXPRESSION "Util Deps: B A SimpleCxx ${simpleCxxDeps} SimpleCxx ${simpleCxxDeps}"
  )

add_test(NAME app_test COMMAND app)
set_tests_properties(app_test PROPERTIES
  PASS_REGULAR_EXPRESSION "Full Deps: WithSubpackages:B A SimpleCxx ${simpleCxxDeps} SimpleCxx ${simpleCxxDeps} A SimpleCxx ${simpleCxxDeps}[;] MixedLang:Mixed Language[;] SimpleCxx:${simpleCxxDeps}"
  )
