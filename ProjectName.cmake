SET(PROJECT_NAME Trilinos)

SET(${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE  ON  CACHE  BOOL
  "Set by default in Trilinos ProjectName.cmake." )

SET(${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE  ON CACHE BOOL
  "Set by default in Trilinos ProjectName.cmake." )

IF ( WIN32 AND NOT CYGWIN )
  MESSAGE(STATUS "Warning: Setting ${PROJECT_NAME}_ENABLE_Fortran=OFF by default"
   " because this is Windows (not cygwin) and we assume to not have Fortran!")
  SET(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT OFF)
ELSE()
  SET(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT ON)
ENDIF()
