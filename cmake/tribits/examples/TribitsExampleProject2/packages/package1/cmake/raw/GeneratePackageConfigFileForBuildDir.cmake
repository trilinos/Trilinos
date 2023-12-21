if (COMMAND tribits_package)
  # Generate Package1Config.cmake file for the build tree (for internal
  # TriBITS-compliant package)
  set(packageBuildDirCMakePackagesDir
    "${${CMAKE_PROJECT_NAME}_BINARY_DIR}/cmake_packages/${PROJECT_NAME}")
  export(EXPORT ${PROJECT_NAME}
    NAMESPACE ${PROJECT_NAME}::
    FILE "${packageBuildDirCMakePackagesDir}/${PROJECT_NAME}ConfigTargets.cmake" )
  configure_file(
    "${CMAKE_CURRENT_LIST_DIR}/Package1Config.cmake.in"
    "${packageBuildDirCMakePackagesDir}/${PROJECT_NAME}/Package1Config.cmake"
    @ONLY )
endif()
