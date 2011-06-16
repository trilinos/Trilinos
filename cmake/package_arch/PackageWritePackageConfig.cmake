INCLUDE(PackageGeneralMacros)

#
#  This function will take a list and turn it into a space separated string
#  adding the prefix to the front of every entry. 
#
FUNCTION(LIST_TO_STRING LIST PREFIX OUTPUT_STRING)
  SET(LIST_STRING "")
  
  FOREACH(ITEM ${LIST})
    SET(LIST_STRING "${LIST_STRING} ${PREFIX}${ITEM}")
  ENDFOREACH()

  SET(${OUTPUT_STRING} ${LIST_STRING} PARENT_SCOPE)
ENDFUNCTION()

#
#  This function will take a list of libraries and turn it into a space
#  separated string. In this case though the prefix is not always added
#  to the front of each entry as libraries can be specified either as a
#  name of a library to find or the absolute path to the library file
#  with any decorations the system uses. When an absolute path is given
#  the entry is used verbatim.
#
FUNCTION(LIBRARY_LIST_TO_STRING LIST PREFIX OUTPUT_STRING)
  SET(LIST_STRING "")
  
  FOREACH(ITEM ${LIST})
    IF(EXISTS ${ITEM})
      SET(LIST_STRING "${LIST_STRING} ${ITEM}")
    ELSE()
      SET(LIST_STRING "${LIST_STRING} ${PREFIX}${ITEM}")
    ENDIF()
  ENDFOREACH()

  SET(${OUTPUT_STRING} ${LIST_STRING} PARENT_SCOPE)
ENDFUNCTION()

#
#  This function checks to see if DEPENDENT_PACKAGE is either a direct or 
#  indirect dependency of PACKAGE_NAME. Optional dependencies are only
#  considered a "dependency" if they are enabled for the package.
#
FUNCTION(CHECK_IS_DEPENDENCY PACKAGE_NAME DEPENDENT_PACKAGE IS_DEPENDENCY)
  #IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
  #  MESSAGE("checking dependency for package ${PACKAGE_NAME} on package ${DEPENDENT_PACKAGE}")
  #ENDIF()
  SET(_IS_DEPENDENCY FALSE)

  #check if the package is being checked against itself
  #a package is always dependent on itself
  IF(PACKAGE_NAME STREQUAL ${DEPENDENT_PACKAGE})
    SET(_IS_DEPENDENCY TRUE)
  ENDIF()

  #check if this is a required dependency
  IF(${PACKAGE_NAME}_LIB_REQUIRED_DEP_PACKAGES)
    LIST(FIND ${PACKAGE_NAME}_LIB_REQUIRED_DEP_PACKAGES ${DEPENDENT_PACKAGE} PACKAGE_IS_REQUIRED_DEP)
    IF(PACKAGE_IS_REQUIRED_DEP GREATER -1)
      SET(_IS_DEPENDENCY TRUE)
    ENDIF()
  ENDIF()

  #check if this is an optional dependency
  IF(${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES)
    LIST(FIND ${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES ${DEPENDENT_PACKAGE} PACKAGE_IS_OPTIONAL_DEP)
    IF(PACKAGE_IS_OPTIONAL_DEP GREATER -1 AND ${PACKAGE_NAME}_ENABLE_${DEPENDENT_PACKAGE})
      SET(_IS_DEPENDENCY TRUE)
    ENDIF()
  ENDIF()

  #if the package is not a direct dependency then test if it is a dependency of the direct
  #dependencies
  IF(NOT _IS_DEPENDENCY)
    #setting to empty because the scoping in cmake has us using the parents list of dependencies
    SET(FULL_DEP_PACKAGES "")
    LIST(APPEND FULL_DEP_PACKAGES ${${PACKAGE_NAME}_LIB_REQUIRED_DEP_PACKAGES} ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES})
    FOREACH(PACKAGE ${FULL_DEP_PACKAGES})
      CHECK_IS_DEPENDENCY(${PACKAGE} ${DEPENDENT_PACKAGE} _IS_INDIRECT_DEPENDENCY)
      IF(_IS_INDIRECT_DEPENDENCY)
        SET(_IS_DEPENDENCY TRUE)
      ENDIF()
    ENDFOREACH()
  ENDIF()

  SET(${IS_DEPENDENCY} ${_IS_DEPENDENCY} PARENT_SCOPE)
ENDFUNCTION()

FUNCTION(PACKAGE_WRITE_PACKAGE_CONFIG_FILE PACKAGE_NAME)

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("For package ${PACKAGE_NAME} creating ${PACKAGE_NAME}Config.cmake")
  ENDIF()

  # Remove from the package list all packages after ${PACKAGE} in the dependency
  # chain. This way each package can create its own minimalist export makefile
  # with no upstream libraries. 
  SET(TMP_PACK_LIST)
  SET(SKIP FALSE)
  FOREACH(PACKAGE ${${PROJECT_NAME}_PACKAGES})
    IF(NOT SKIP)
      LIST(APPEND TMP_PACK_LIST ${PACKAGE})
      IF(PACKAGE_NAME STREQUAL ${PACKAGE})
        SET(SKIP TRUE)
      ENDIF()
    ENDIF()
  ENDFOREACH()

  SET(PACKAGE_LIST ${TMP_PACK_LIST})

  # Reverse the order of the package list, letting us loop 
  # from most-dependent to least-dependent. 
  LIST(REVERSE PACKAGE_LIST)

  # Now that we have a reduced set of packages check each one to see if A. PACKAGE_NAME actually
  # depends on PACKAGE and B. PACKAGE is enabled (since it could be an optional dependency.) If
  # both A and B are true then we add their libraries to the list and to the list of packages so
  # we can then loop over their tpls later

  SET(FULL_PACKAGE_SET "")
  SET(FULL_LIBRARY_SET "")
  FOREACH(PACKAGE ${PACKAGE_LIST})
    CHECK_IS_DEPENDENCY(${PACKAGE_NAME} ${PACKAGE} IS_DEPENDENCY)
    
    IF(IS_DEPENDENCY)
      LIST(APPEND FULL_PACKAGE_SET ${PACKAGE})
      LIST(APPEND FULL_LIBRARY_SET ${${PACKAGE}_LIBRARIES})
    ENDIF()
  ENDFOREACH()

  SET(FULL_TPL_SET "")
  FOREACH(PACKAGE ${FULL_PACKAGE_SET})
    LIST(APPEND FULL_TPL_SET ${${PACKAGE}_LIB_REQUIRED_DEP_TPLS})

    SET(OPTIONAL_TPLS ${${PACKAGE}_LIB_OPTIONAL_DEP_TPLS})

    FOREACH(TPL ${OPTIONAL_TPLS})
      IF(${PACKAGE}_ENABLE_${TPL})
        LIST(APPEND FULL_TPL_SET ${TPL})
      ENDIF()
    ENDFOREACH()
  ENDFOREACH()

  #We will use the complete list of supported tpls for the project
  #to help us create a properly ordered list of tpls.
  IF(FULL_TPL_SET)
    # Reversing the tpl list so that the list of tpls will be produced in
    # order of most dependent to least dependent.
    SET(TPL_LIST ${${PROJECT_NAME}_TPLS})
    LIST(REVERSE TPL_LIST)
    
    SET(ORDERED_FULL_TPL_SET "")
    FOREACH(TPL ${TPL_LIST})
      LIST(FIND FULL_TPL_SET ${TPL} FOUND)
      
      IF(FOUND GREATER -1)
        LIST(APPEND ORDERED_FULL_TPL_SET ${TPL})
      ENDIF()
    ENDFOREACH()
  ENDIF()
  
  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("ORDERED_FULL_TPL_SET after all tpls = ${ORDERED_FULL_TPL_SET}")
  ENDIF()
  
  SET(${PACKAGE_NAME}_TPL_LIBRARIES "")
  SET(${PACKAGE_NAME}_TPL_INCLUDE_DIRS "")
  SET(${PACKAGE_NAME}_TPL_LIBRARY_DIRS "")
  FOREACH(TPL ${ORDERED_FULL_TPL_SET})
    LIST(APPEND ${PACKAGE_NAME}_TPL_LIBRARIES ${TPL_${TPL}_LIBRARIES})
    LIST(APPEND ${PACKAGE_NAME}_TPL_INCLUDE_DIRS ${TPL_${TPL}_INCLUDE_DIRS})
    LIST(APPEND ${PACKAGE_NAME}_TPL_LIBRARY_DIRS ${TPL_${TPL}_LIBRARY_DIRS})
  ENDFOREACH()
  
  # Generate a note discouraging editing of the <package>Config.cmake file
  SET(DISCOURAGE_EDITING "Do not edit: This file was generated automatically by CMake.")
  
  ######
  # Create a configure file for the build tree. Creating this file in the base dir
  # of the package since it is possible that the cmake returned path where the file
  # was found would be useful for a package, and having to dig through the hiding that
  # is done wouldn't be nice.
  ######

  SET(LIBRARY_DIRS ${${PACKAGE_NAME}_LIBRARY_DIRS})
  SET(INCLUDE_DIRS ${${PACKAGE_NAME}_INCLUDE_DIRS})

  # Custom code in configuration file.
  SET(PACKAGE_CONFIG_CODE "")

  # Import build tree targets into applications.
  IF(FULL_LIBRARY_SET)
    SET(PACKAGE_CONFIG_CODE "${PACKAGE_CONFIG_CODE}
# Import ${PROJECT_NAME} targets
IF(NOT ${PROJECT_NAME}_TARGETS_IMPORTED)
  SET(${PROJECT_NAME}_TARGETS_IMPORTED 1)
  INCLUDE(\"${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake\")
ENDIF()
")
  ENDIF()

  # Write the specification of the rpath if necessary. This is only needed if we're building shared libraries. 
  IF(BUILD_SHARED_LIBS)
    STRING(REPLACE ";" ":" SHARED_LIB_RPATH_COMMAND "${LIBRARY_DIRS}")
    SET(SHARED_LIB_RPATH_COMMAND ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${SHARED_LIB_RPATH_COMMAND})
  ENDIF()

  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/cmake/PackageConfig.cmake.in 
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}Config.cmake
    )

  IF(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)
    ######
    # Create a Makefile.export.<package_name> for the build tree. This is the equivalent
    # of the cmake version only slightly changed so that it can be directly imported into
    # a Makefile.
    ######

    LIST_TO_STRING("${FULL_LIBRARY_SET}" ${CMAKE_LINK_LIBRARY_FLAG} MAKEFILE_FULL_LIBRARY_SET)
    LIST_TO_STRING("${LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_LIBRARY_DIRS)
    LIST_TO_STRING("${INCLUDE_DIRS}" "-I" MAKEFILE_INCLUDE_DIRS)
    LIST_TO_STRING("${${PACKAGE_NAME}_TPL_INCLUDE_DIRS}" "-I" MAKEFILE_${PACKAGE_NAME}_TPL_INCLUDE_DIRS)
    LIST_TO_STRING("${${PACKAGE_NAME}_TPL_LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_${PACKAGE_NAME}_TPL_LIBRARY_DIRS)
    #the TPL library names have to be treated differently
    LIBRARY_LIST_TO_STRING("${${PACKAGE_NAME}_TPL_LIBRARIES}" ${CMAKE_LINK_LIBRARY_FLAG} MAKEFILE_${PACKAGE_NAME}_TPL_LIBRARIES)

    LIBRARY_LIST_TO_STRING("${${TPL_MPI_LIBRARIES}}" ${CMAKE_LINK_LIBRARY_FLAG} "MAKEFILE_TPL_MPI_LIBRARIES")
    LIST_TO_STRING("${${TPL_MPI_LIBRARY_DIRS}}" ${CMAKE_LIBRARY_PATH_FLAG} "MAKEFILE_TPL_MPI_LIBRARY_DIRS")
    LIST_TO_STRING("${${TPL_MPI_INCLUDE_DIRS}}" "-I" "MAKEFILE_TPL_MPI_INCLUDE_DIRS")
    
    LIST_TO_STRING("${FULL_PACKAGE_SET}" "" MAKEFILE_FULL_PACKAGE_SET)
    LIST_TO_STRING("${ORDERED_FULL_TPL_SET}" "" MAKEFILE_ORDERED_FULL_TPL_SET)

    #create an upper case name of the package so that we can make deprecated versions of them to help people
    #transistioning from the autotools version diagnose any missed variables.
    STRING(TOUPPER ${PACKAGE_NAME} PACKAGE_NAME_UPPER)

    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/cmake/PackageConfig.export.in 
      ${CMAKE_CURRENT_BINARY_DIR}/Makefile.export.${PACKAGE_NAME}
      )
  ENDIF()

  ######
  # Create a configure file for the install tree and set the install target for it. This
  # file isn't generally useful inside the build tree so it is being "hidden" in the 
  # CMakeFiles directory. It will be placed in the base install directory for ${PROJECT_NAME}
  # when installed.
  ######

  SET(LIBRARY_DIRS ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})
  SET(INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

  # Custom code in configuration file.
  SET(PACKAGE_CONFIG_CODE "")

  # Import install tree targets into applications.
  GET_PROPERTY(HAS_INSTALL_TARGETS GLOBAL PROPERTY ${PACKAGE_NAME}_HAS_INSTALL_TARGETS)
  IF(HAS_INSTALL_TARGETS)
    SET(PACKAGE_CONFIG_CODE "${PACKAGE_CONFIG_CODE}
# Import ${PROJECT_NAME} targets
IF(NOT ${PROJECT_NAME}_TARGETS_IMPORTED)
  SET(${PROJECT_NAME}_TARGETS_IMPORTED 1)
  INCLUDE(\"${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PROJECT_NAME}/${PROJECT_NAME}Targets.cmake\")
ENDIF()
")
  ENDIF()

  # Write the specification of the rpath if necessary. This is only needed if we're building shared libraries. 
  IF(BUILD_SHARED_LIBS)
    SET(SHARED_LIB_RPATH_COMMAND ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})
  ENDIF()

  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/cmake/PackageConfig.cmake.in 
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${PACKAGE_NAME}Config_install.cmake
    )

  INSTALL(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${PACKAGE_NAME}Config_install.cmake
    DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PACKAGE_NAME}"
    RENAME ${PACKAGE_NAME}Config.cmake
  )

  IF(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)
    ######
    # Create a Makefile.export.<package_name> for the install tree. This is the equivalent
    # of the cmake version only slightly changed so that it can be directly imported into
    # a Makefile.
    ######

    LIST_TO_STRING("${LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_LIBRARY_DIRS)
    LIST_TO_STRING("${INCLUDE_DIRS}" "-I" MAKEFILE_INCLUDE_DIRS)

    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/cmake/PackageConfig.export.in 
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Makefile.export.${PACKAGE_NAME}_install
      )

    INSTALL(
      FILES ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Makefile.export.${PACKAGE_NAME}_install
      DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}"
      RENAME Makefile.export.${PACKAGE_NAME}
    )
  ENDIF()
ENDFUNCTION()

FUNCTION(PACKAGE_ARCH_WRITE_CONFIG_FILE)

  # Reversing the package list so that libraries will be produced in order of
  # most dependent to least dependent.
  SET(PACKAGE_LIST ${${PROJECT_NAME}_PACKAGES})
  LIST(REVERSE PACKAGE_LIST)
  

  # Loop over all packages to determine which were enabled. Then build a list
  # of all their libraries/includes in the proper order for linking
  SET(FULL_PACKAGE_SET "")
  SET(FULL_LIBRARY_SET "")
  SET(FULL_INCLUDE_DIRS_SET "")
  SET(FULL_LIBRARY_DIRS_SET "")
  FOREACH(PACKAGE ${PACKAGE_LIST})
    IF(${PROJECT_NAME}_ENABLE_${PACKAGE})
      LIST(APPEND FULL_PACKAGE_SET ${PACKAGE})
      LIST(APPEND FULL_LIBRARY_SET ${${PACKAGE}_LIBRARIES})
      LIST(APPEND FULL_INCLUDE_DIRS_SET ${${PACKAGE}_INCLUDE_DIRS})
      LIST(APPEND FULL_LIBRARY_DIRS_SET ${${PACKAGE}_LIBRARY_DIRS})
    ENDIF()
  ENDFOREACH()
  
  SET(${PROJECT_NAME}_CONFIG_LIBRARIES ${FULL_LIBRARY_SET})
  
  # Reversing the tpl list so that the list of tpls will be produced in
  # order of most dependent to least dependent.
  SET(TPL_LIST ${${PROJECT_NAME}_TPLS})
  LIST(REVERSE TPL_LIST)
  
  # Loop over all TPLs to determine which were enabled. Then build a list
  # of all their libraries/includes in the proper order for linking
  SET(FULL_TPL_SET "")
  SET(FULL_TPL_LIBRARY_SET "")
  SET(FULL_TPL_INCLUDE_DIRS_SET "")
  SET(FULL_TPL_LIBRARY_DIRS_SET "")
  FOREACH(TPL ${TPL_LIST})
    IF(TPL_ENABLE_${TPL})
      LIST(APPEND FULL_TPL_SET ${TPL})
      LIST(APPEND FULL_TPL_LIBRARY_SET ${TPL_${TPL}_LIBRARIES})
      LIST(APPEND FULL_TPL_INCLUDE_DIRS_SET ${TPL_${TPL}_INCLUDE_DIRS})
      LIST(APPEND FULL_TPL_LIBRARY_DIRS_SET ${TPL_${TPL}_LIBRARY_DIRS})
    ENDIF()
  ENDFOREACH()
  
  # it is possible that tpls are in the same directory, to keep from 
  # having a very long include path or library path we will strip out
  # any duplicates. This shouldn't affect which include or library is
  # found since the first instance of any path will be the one that is
  # kept.
  LIST(REMOVE_DUPLICATES FULL_TPL_INCLUDE_DIRS_SET)
  LIST(REMOVE_DUPLICATES FULL_TPL_LIBRARY_DIRS_SET)
  
  SET(${PROJECT_NAME}_CONFIG_TPL_INCLUDE_DIRS ${FULL_TPL_INCLUDE_DIRS_SET})
  SET(${PROJECT_NAME}_CONFIG_TPL_LIBRARY_DIRS ${FULL_TPL_LIBRARY_DIRS_SET})
  SET(${PROJECT_NAME}_CONFIG_TPL_LIBRARIES ${FULL_TPL_LIBRARY_SET})
  
  #
  # Configure two files for finding ${PROJECT_NAME}. One for the build tree and one for installing
  #

  # Generate a note discouraging editing of the <package>Config.cmake file
  SET(DISCOURAGE_EDITING "Do not edit: This file was generated automatically by CMake.")
  
  #Config file for setting variables and finding include/library paths from the build directory
  SET(${PROJECT_NAME}_CONFIG_INCLUDE_DIRS ${FULL_INCLUDE_DIRS_SET})
  SET(${PROJECT_NAME}_CONFIG_LIBRARY_DIRS ${FULL_LIBRARY_DIRS_SET})

  # Write the specification of the rpath if necessary. This is only needed if we're building shared libraries. 
  IF(BUILD_SHARED_LIBS)
    STRING(REPLACE ";" ":" SHARED_LIB_RPATH_COMMAND "${${PROJECT_NAME}_CONFIG_LIBRARY_DIRS}")
    SET(SHARED_LIB_RPATH_COMMAND ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${SHARED_LIB_RPATH_COMMAND})
  ENDIF()

  # Custom code in configuration file.
  SET(PROJECT_CONFIG_CODE "")

  # Export targets from the build tree.
  IF(FULL_LIBRARY_SET)
    EXPORT(TARGETS ${FULL_LIBRARY_SET} FILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake")
    # Import the targets in applications.
    SET(PROJECT_CONFIG_CODE "${PROJECT_CONFIG_CODE}
# Import ${PROJECT_NAME} targets
IF(NOT ${PROJECT_NAME}_TARGETS_IMPORTED)
  SET(${PROJECT_NAME}_TARGETS_IMPORTED 1)
  INCLUDE(\"${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake\")
ENDIF()
")
  ENDIF()

  # Appending the logic to include each package's config file.
  SET(LOAD_CODE "# Load configurations from enabled packages\n")
  FOREACH(PACKAGE ${FULL_PACKAGE_SET})
    SET(LOAD_CODE "${LOAD_CODE}include(\"${${PACKAGE}_BINARY_DIR}/${PACKAGE}Config.cmake\")\n")
  ENDFOREACH()
  SET(PROJECT_CONFIG_CODE "${PROJECT_CONFIG_CODE}\n${LOAD_CODE}")

  CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake )

  IF(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)
    ######
    # Create a Makefile.export.<project_name> for the build tree. This is the equivalent
    # of the cmake version only slightly changed so that it can be directly imported into
    # a Makefile.
    ######

    LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_LIBRARIES}" ${CMAKE_LINK_LIBRARY_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_LIBRARIES)
    LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_LIBRARY_DIRS)
    LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_INCLUDE_DIRS}" "-I" MAKEFILE_${PROJECT_NAME}_CONFIG_INCLUDE_DIRS)
    LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_TPL_INCLUDE_DIRS}" "-I" MAKEFILE_${PROJECT_NAME}_CONFIG_TPL_INCLUDE_DIRS)
    LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_TPL_LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_TPL_LIBRARY_DIRS)
    #the TPL library names have to be treated differently
    LIBRARY_LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_TPL_LIBRARIES}" ${CMAKE_LINK_LIBRARY_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_TPL_LIBRARIES)

    LIBRARY_LIST_TO_STRING("${${TPL_MPI_LIBRARIES}}" ${CMAKE_LINK_LIBRARY_FLAG} "MAKEFILE_TPL_MPI_LIBRARIES")
    LIST_TO_STRING("${${TPL_MPI_LIBRARY_DIRS}}" ${CMAKE_LIBRARY_PATH_FLAG} "MAKEFILE_TPL_MPI_LIBRARY_DIRS")
    LIST_TO_STRING("${${TPL_MPI_INCLUDE_DIRS}}" "-I" "MAKEFILE_TPL_MPI_INCLUDE_DIRS")
    
    LIST_TO_STRING("${FULL_PACKAGE_SET}" "" MAKEFILE_FULL_PACKAGE_SET)
    LIST_TO_STRING("${FULL_TPL_SET}" "" MAKEFILE_FULL_TPL_SET)

    CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.export.in
      ${CMAKE_CURRENT_BINARY_DIR}/Makefile.export.${PROJECT_NAME})
  ENDIF()  

  ######
  # Create a configure file for the install tree and set the install target for it. This
  # file isn't generally useful inside the build tree. It will be placed in the base
  # install directory for ${PROJECT_NAME} when installed.
  ######

  SET(${PROJECT_NAME}_CONFIG_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
  SET(${PROJECT_NAME}_CONFIG_LIBRARY_DIRS ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})

  # Write the specification of the rpath if necessary. This is only needed if we're building shared libraries. 
  IF(BUILD_SHARED_LIBS)
    SET(SHARED_LIB_RPATH_COMMAND ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})
  ENDIF()

  # Custom code in configuration file.
  SET(PROJECT_CONFIG_CODE "")

  # Export targets from the install tree.
  GET_PROPERTY(HAS_INSTALL_TARGETS GLOBAL PROPERTY ${PROJECT_NAME}_HAS_INSTALL_TARGETS)
  IF(HAS_INSTALL_TARGETS)
    INSTALL(
      EXPORT ${PROJECT_NAME}
      DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PROJECT_NAME}"
      FILE ${PROJECT_NAME}Targets.cmake
      )
    # Import the targets in applications.
    SET(PROJECT_CONFIG_CODE "${PROJECT_CONFIG_CODE}
# Import ${PROJECT_NAME} targets
IF(NOT ${PROJECT_NAME}_TARGETS_IMPORTED)
  SET(${PROJECT_NAME}_TARGETS_IMPORTED 1)
  INCLUDE(\"${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PROJECT_NAME}/${PROJECT_NAME}Targets.cmake\")
ENDIF()
")
  ENDIF()

  # Appending the logic to include each package's config file.
  SET(LOAD_CODE "# Load configurations from enabled packages\n")
  FOREACH(PACKAGE ${FULL_PACKAGE_SET})
    SET(LOAD_CODE "${LOAD_CODE}include(\"${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PACKAGE}/${PACKAGE}Config.cmake\")\n")
  ENDFOREACH()
  SET(PROJECT_CONFIG_CODE "${PROJECT_CONFIG_CODE}\n${LOAD_CODE}")

  CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake )

  INSTALL(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake
    DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PROJECT_NAME}"
    RENAME ${PROJECT_NAME}Config.cmake
  )
  
  IF(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)
    ######
    # Create a Makefile.export.<project_name> for the install tree. This is the equivalent
    # of the cmake version only slightly changed so that it can be directly imported into
    # a Makefile.
    ######

    LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_LIBRARY_DIRS)
    LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_INCLUDE_DIRS}" "-I" MAKEFILE_${PROJECT_NAME}_CONFIG_INCLUDE_DIRS)

    CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.export.in
      ${CMAKE_CURRENT_BINARY_DIR}/Makefile.export.${PROJECT_NAME}_install )

    INSTALL(
      FILES ${CMAKE_CURRENT_BINARY_DIR}/Makefile.export.${PROJECT_NAME}_install
      DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}"
      RENAME Makefile.export.${PROJECT_NAME}
    )
  ENDIF()
  
  #
  # Configure the version file for ${PROJECT_NAME}
  #
  
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}ConfigVersion.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake @ONLY)

  INSTALL(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PROJECT_NAME}"
  )

ENDFUNCTION()

