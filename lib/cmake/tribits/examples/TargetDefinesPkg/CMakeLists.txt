tribits_package(TargetDefinesPkg)

tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

tribits_add_library(
  targetdefinespkg
  HEADERS TargetDefinesPkg.hpp
  SOURCES TargetDefinesPkg.cpp
  )

tribits_add_test_directories(test)

tribits_package_postprocess()
