# This file is meant to be included in the configuration of Trilinos to
# disable packages internally.  These excludes can be overwritten on the
# command line in the CMake cache.

INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/casl-exclude-trilinos-packages.cmake)

FOREACH(EXCLUDED_PACKAGE ${Trilinos_EXCLUDE_PACKAGES})
  SET(${PROJECT_NAME}_ENABLE_${EXCLUDED_PACKAGE} OFF CACHE BOOL
    "Disabled in casl-disable-excluded-trilinos-packages.cmake")
ENDFOREACH()
