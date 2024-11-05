# Generate and install the Package1Config.cmake file for the install tree
# (needed for both internal and external TriBITS package)
set(pkgConfigInstallDir "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
install(EXPORT ${PROJECT_NAME}
  DESTINATION "${pkgConfigInstallDir}"
  NAMESPACE ${PROJECT_NAME}::
  FILE ${PROJECT_NAME}ConfigTargets.cmake )
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/raw/Package1Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Package1Config.install.cmake"
  @ONLY )
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Package1Config.install.cmake"
  RENAME "Package1Config.cmake"
  DESTINATION "${pkgConfigInstallDir}" )
