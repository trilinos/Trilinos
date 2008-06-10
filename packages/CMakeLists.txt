# $Header$

# 2008/06/09: rabartl: The list of packages should probably be defined
#   in this file and not in the file above this one.

FOREACH(PACKAGE ${${PROJECT_NAME}_PACKAGES})
  IF(${${PROJECT_NAME}_ENABLE_${PACKAGE}} STREQUAL "ON")
    STRING(TOLOWER ${PACKAGE} PACKAGE_DIR)
    # 2008/06/09: rabartl: Above, you can not assume that all package names are lower
    #   case.  This will not work for PyTrilinos, ForTrilinos, ThreadPool, etc.  You should
    #   just name the packages with their directory names!
    ADD_SUBDIRECTORY(${PACKAGE_DIR})
  ENDIF(${${PROJECT_NAME}_ENABLE_${PACKAGE}} STREQUAL "ON")
ENDFOREACH(PACKAGE ${PROJECT_NAME}_PACKAGES)
