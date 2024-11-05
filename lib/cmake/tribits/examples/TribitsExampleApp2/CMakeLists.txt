cmake_minimum_required(VERSION 3.23.1)

project(TribitsExApp2
  DESCRIPTION
    "Example raw CMake project using packages installed from TribitsExampleProject2"
  VERSION 0.0.0
  LANGUAGES NONE  # Defined below after reading in compilers
  )

set(${PROJECT_NAME}_USE_COMPONENTS "" CACHE STRING
  "Components/Packages to use from TribitsExampleProject2: <C0>,<C1>,..." )
string(REPLACE "," ";" ${PROJECT_NAME}_USE_COMPONENTS
  "${${PROJECT_NAME}_USE_COMPONENTS}")

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
include(AppHelperFuncs)

getTribitsExProj2StuffForApp()

# Enable the compilers now that we have gotten them from the *Config.cmake file
enable_language(C)
enable_language(CXX)
if (CMAKE_Fortran_COMPILER)
  enable_language(Fortran)
endif()

# Build the APP and link to libraries from TribitsExProj packages
add_executable(app app.cpp)
target_link_libraries(app PRIVATE ${APP_DEPS_LIB_TARGETS})
addAppDepCompileDefines()

# Set up tests

enable_testing()

getExpectedAppDepsStr(expectedDepsStr)

add_test(app_test app)
set_tests_properties(app_test PROPERTIES
  PASS_REGULAR_EXPRESSION "Full Deps: ${expectedDepsStr}"
  )
