cmake_minimum_required(VERSION 3.17.0)

add_library(simpletpl SimpleTpl.hpp SimpleTpl.cpp)

install(
  TARGETS simpletpl
  RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
  )

install(
  FILES SimpleTpl.hpp
  DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
  )
