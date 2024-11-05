cmake_minimum_required(VERSION 3.23.1)

project(TribitsOldSimpleExApp
  DESCRIPTION
    "Example raw CMake project using packages installed from new or old TribitsExampleProject"
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
if (${PROJECT_NAME}_USE_DEPRECATED_TARGETS)
  target_link_libraries(util PRIVATE pws_b)
  # NOTE: Above will generate a DEPRECATION warning for newer TriBITS
  # versions!
else()
  target_link_libraries(util PRIVATE ${WithSubpackagesB_LIBRARIES})
  # NOTE: Above is the correct way to link against libraries from an
  # individual package for old and new TriBITS!
endif()
target_include_directories(util
  PRIVATE ${SimpleCxx_INCLUDE_DIRS} ${SimpleCxx_TPL_INCLUDE_DIRS})

# Build the APP and link to libraries from TribitsExProj packages
add_executable(app app.cpp)
target_link_libraries(app
  PRIVATE ${TribitsExProj_LIBRARIES} ${TribitsExProj_TPL_LIBRARIES})
target_include_directories(app
  PRIVATE ${TribitsExProj_INCLUDE_DIRS} ${TribitsExProj_TPL_INCLUDE_DIRS})

# Set up tests

enable_testing()

if ("SimpleTpl"  IN_LIST  TribitsExProj_TPL_LIST  OR  TARGET  SimpleTpl::all_libs)
  set(simpleCxxDeps "simpletpl headeronlytpl")
else()
  set(simpleCxxDeps "headeronlytpl")
endif()
# NOTE: Above, TribitsExProj_TPL_LIST does not exist in new TriBITS so to make
# this work for new TriBITS as well, we need to check for target
# SimpleTpl::all_libs.  Small price to pay for progress!

add_test(NAME util_test COMMAND util)
set_tests_properties(util_test PROPERTIES
  PASS_REGULAR_EXPRESSION "Util Deps: B A SimpleCxx ${simpleCxxDeps} SimpleCxx ${simpleCxxDeps}"
  )

add_test(NAME app_test COMMAND app)
set_tests_properties(app_test PROPERTIES
  PASS_REGULAR_EXPRESSION "Full Deps: WithSubpackages:B A SimpleCxx ${simpleCxxDeps} SimpleCxx ${simpleCxxDeps} A SimpleCxx ${simpleCxxDeps}[;] MixedLang:Mixed Language[;] SimpleCxx:${simpleCxxDeps}"
  )
