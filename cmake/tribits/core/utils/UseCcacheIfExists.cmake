# Include this file in your CMake project with either -C
# UseCcacheIfExists.cmake or with the TriBITS option -D
# <Project>_CONFIGURE_OPTIONS_FILE=UseCcacheIfExists.cmake.

FIND_PROGRAM(CCACHE_EXEC ccache)
IF (CCACHE_EXEC)
  MESSAGE("-- NOTE: ccache program exists so using it for RULE_LAUNCH_COMPILE!") 
  SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_EXEC}")
ENDIF()
