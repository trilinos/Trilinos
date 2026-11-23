# Ninja build generator configuration

# Use Ninja instead of Unix Makefiles for faster parallel builds
set(CMAKE_GENERATOR "Ninja" CACHE STRING "")

# Colored diagnostic output (works in most terminals)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always" CACHE STRING "" FORCE)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdiagnostics-color=always" CACHE STRING "" FORCE)
