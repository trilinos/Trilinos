cmake_minimum_required(VERSION 2.7.20090924)

# Include the common body provided by TriBITS. This should be all that
# is necessary.
set(TRIBITS_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../tribits)
include(${TRIBITS_ROOT}/ctest/tdd/TribitsDriverCMakeLists.cmake)