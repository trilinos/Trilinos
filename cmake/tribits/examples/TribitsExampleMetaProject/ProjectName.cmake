# Must set the project name at very beginning before including anything else
SET(PROJECT_NAME TribitsExMetaProj)

# Turn on export dependency generation for WrapExteranl package
SET(${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES_DEFAULT ON)

# Must always have the PRE extra repo TribitsExampleProject/ or can't build
SET(${PROJECT_NAME}_EXTRAREPOS_FILE
  ${CMAKE_CURRENT_LIST_DIR}/cmake/ExtraRepositoriesList.cmake
  CACHE  FILEPATH  "Set in ProjectName.cmake")
SET(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE  Continuous
  CACHE  STRING  "Set in ProjectName.cmake")

# Generate <Project>RepoVersion.txt file
SET(${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE_DEFAULT ON)

# Do the all-at-once approach with tribits_ctest_driver() by default
SET(${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE_DEFAULT ON)
