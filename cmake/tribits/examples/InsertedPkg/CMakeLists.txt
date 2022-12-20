tribits_package(InsertedPkg)

tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

tribits_add_library(
  externalpkg
  HEADERS InsertedPkg.hpp
  SOURCES InsertedPkg.cpp
  )

tribits_add_test_directories(test)

tribits_package_postprocess()
