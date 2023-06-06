include_raw_cmake_build()

tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

tribits_add_library(package1  HEADERS  Package1.hpp  SOURCES  Package1.cpp)

tribits_add_executable(package1-prg  NOEXEPREFIX  NOEXESUFFIX
  SOURCES  Package1_Prg.cpp  INSTALLABLE )
