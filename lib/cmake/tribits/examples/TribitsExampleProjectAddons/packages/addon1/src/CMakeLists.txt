tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})
tribits_add_library(addon1
  SOURCES Addon1.cpp
  HEADERS Addon1.hpp
  NOINSTALLHEADERS 
)
