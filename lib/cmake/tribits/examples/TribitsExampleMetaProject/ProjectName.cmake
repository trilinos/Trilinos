# Must set the project name at very beginning before including anything else
set(PROJECT_NAME TribitsExMetaProj)

# Turn on export dependency generation for WrapExteranl package
set(${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES_DEFAULT ON)

# Must always have the PRE extra repo TribitsExampleProject/ or can't build
set(${PROJECT_NAME}_EXTRAREPOS_FILE
  ${CMAKE_CURRENT_LIST_DIR}/cmake/ExtraRepositoriesList.cmake
  CACHE  FILEPATH  "Set in ProjectName.cmake")
set(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE  Continuous
  CACHE  STRING  "Set in ProjectName.cmake")

# Generate <Project>RepoVersion.txt file
set(${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE_DEFAULT ON)

# Do the all-at-once approach with tribits_ctest_driver() by default
set(${PROJECT_NAME}_CTEST_DO_ALL_AT_ONCE_DEFAULT ON)
