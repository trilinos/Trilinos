# Generate the all_libs target(s)
add_library(Package1_all_libs INTERFACE)
set_target_properties(Package1_all_libs
  PROPERTIES EXPORT_NAME all_libs)
target_link_libraries(Package1_all_libs
  INTERFACE Package1_package1)
install(TARGETS Package1_all_libs
  EXPORT ${PROJECT_NAME}
  COMPONENT ${PROJECT_NAME}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )
add_library(Package1::all_libs ALIAS Package1_all_libs)
