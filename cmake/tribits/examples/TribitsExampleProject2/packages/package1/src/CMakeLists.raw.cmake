# Create and install library 'package1'
add_library(Package1_package1 Package1.hpp Package1.cpp)
target_include_directories(Package1_package1
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(Package1_package1
  PRIVATE tpl1::tpl1 )
set_target_properties(Package1_package1 PROPERTIES
  EXPORT_NAME package1)
add_library(Package1::package1 ALIAS Package1_package1)
install(
  TARGETS Package1_package1
  EXPORT ${PROJECT_NAME}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )
install(
  FILES Package1.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )

# Create and install executable 'package1-prg' 
add_executable(package1-prg Package1_Prg.cpp)
target_link_libraries(package1-prg PRIVATE Package1::package1)
install(
  TARGETS package1-prg
  EXPORT ${PROJECT_NAME}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )
