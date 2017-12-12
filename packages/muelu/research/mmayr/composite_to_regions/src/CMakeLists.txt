INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_GLOB(SOURCES ${DIR}/*.cpp)

#
# Executables
#

  TRIBITS_ADD_EXECUTABLE(
    composite_to_region_driver 
    SOURCES main.cpp
  )

