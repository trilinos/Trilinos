# Include this file in your CMake project with either -C
# UseCcacheIfExists.cmake or with the TriBITS option -D
# <Project>_CONFIGURE_OPTIONS_FILE=UseCcacheIfExists.cmake.

find_program(CCACHE_EXEC ccache)
if (CCACHE_EXEC)
  message("-- NOTE: ccache program exists so using it for RULE_LAUNCH_COMPILE!") 
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_EXEC}")
endif()
