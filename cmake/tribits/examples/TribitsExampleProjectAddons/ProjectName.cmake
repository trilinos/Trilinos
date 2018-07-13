# Must set the project name at very beginning before including anything else
SET(PROJECT_NAME TribitsExProjAddons)

# Turn on export dependency generation for WrapExteranl package
SET(${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES_DEFAULT ON)

# Must always have the PRE extra repo TribitsExampleProject/ or can't build
SET(${PROJECT_NAME}_EXTRAREPOS_FILE  cmake/ExtraRepositoriesList.cmake
  CACHE  FILEPATH  "Set in ProjectName.cmake")
SET(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE  Continuous
  CACHE  STRING  "Set in ProjectName.cmake")
