INCLUDE(PackageGeneralMacros)

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

  SET(FULL_PACKAGE_SET ${PACKAGE_NAME})
  SET(FULL_LIBRARY_SET "")
  LIST(APPEND FULL_LIBRARY_SET ${${PACKAGE_NAME}_LIBRARIES})
  FOREACH(PACKAGE ${PACKAGE_LIST})
    IF(${PACKAGE_NAME}_LIB_REQUIRED_DEP_PACKAGES)
      LIST(FIND ${PACKAGE_NAME}_LIB_REQUIRED_DEP_PACKAGES ${PACKAGE} PACKAGE_IS_REQUIRED_DEP)
    ENDIF()
    IF(${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES)
      LIST(FIND ${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES ${PACKAGE} PACKAGE_IS_OPTIONAL_DEP)
    ENDIF()

    IF(PACKAGE_IS_REQUIRED_DEP GREATER -1)
      LIST(APPEND FULL_PACKAGE_SET ${PACKAGE})
      LIST(APPEND FULL_LIBRARY_SET ${${PACKAGE}_LIBRARIES})
    ENDIF()
    IF(PACKAGE_IS_OPTIONAL_DEP GREATER -1 AND ${PACKAGE_NAME}_ENABLE_${PACKAGE})
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
  
  # Create a configure file for the build tree. Creating this file in the base dir
  # of the package since it is possible that the cmake returned path where the file
  # was found would be useful for a package, and having to dig through the hiding that
  # is done wouldn't be nice.
  SET(LIBRARY_DIRS ${${PACKAGE_NAME}_LIBRARY_DIRS})
  SET(INCLUDE_DIRS ${${PACKAGE_NAME}_INCLUDE_DIRS})

  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/cmake/PackageConfig.cmake.in 
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}Config.cmake
    )

  # Create a configure file for the install tree and set the install target for it. This
  # file isn't generally useful inside the build tree so it is being "hidden" in the 
  # CMakeFiles directory. It will be placed in the base install directory for Trilinos
  # when installed.
  SET(LIBRARY_DIRS ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})
  SET(INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

  CONFIGURE_FILE(
    ${PROJECT_SOURCE_DIR}/cmake/PackageConfig.cmake.in 
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${PACKAGE_NAME}Config_install.cmake
    )

  INSTALL(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${PACKAGE_NAME}Config_install.cmake
    DESTINATION "."
    RENAME ${PACKAGE_NAME}Config.cmake
  )

ENDFUNCTION()

FUNCTION(PACKAGE_WRITE_TRILINOS_CONFIG_FILE)

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
  # Configure two files for finding Trilinos. One for the build tree and one for installing
  #

  # Generate a note discouraging editing of the <package>Config.cmake file
  SET(DISCOURAGE_EDITING "Do not edit: This file was generated automatically by CMake.")
  
 #Config file for setting variables and finding include/library paths from the build directory
  SET(${PROJECT_NAME}_CONFIG_INCLUDE_DIRS ${FULL_INCLUDE_DIRS_SET})
  SET(${PROJECT_NAME}_CONFIG_LIBRARY_DIRS ${FULL_LIBRARY_DIRS_SET})

  CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/TrilinosConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/TrilinosConfig.cmake )

  # To be able to properly pull in the configure files for each package we have to append some cmake code
  # to the end of the configured file. This is that code. The "FOREACH_BLOCK" is same for both the install
  # and build tree configure files, however, the files that are globbed are different. 
  SET(GLOB_LINE "FILE(GLOB PACKAGE_CONFIG_FILES \"${CMAKE_CURRENT_BINARY_DIR}/packages/*/*Config.cmake\")\n")
  SET(FOREACH_BLOCK "FOREACH(FILE \${PACKAGE_CONFIG_FILES})\n  IF(NOT \${FILE} MATCHES \"${PROJECT_NAME}Config.cmake\")\n    INCLUDE(\${FILE})\n  ENDIF()\nENDFOREACH()\n")
  FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/TrilinosConfig.cmake ${GLOB_LINE})
  FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/TrilinosConfig.cmake ${FOREACH_BLOCK})

  #Config file for setting variables and finding include/library paths from the Install directory
  SET(${PROJECT_NAME}_CONFIG_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
  SET(${PROJECT_NAME}_CONFIG_LIBRARY_DIRS ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})

  CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/TrilinosConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/TrilinosConfig_install.cmake )

  # Appending the logic to include each package's config file.
  SET(GLOB_LINE "FILE(GLOB PACKAGE_CONFIG_FILES \"${CMAKE_INSTALL_PREFIX}/*Config.cmake\")\n")
  FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/TrilinosConfig_install.cmake ${GLOB_LINE})
  FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/TrilinosConfig_install.cmake ${FOREACH_BLOCK})

  INSTALL(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/TrilinosConfig_install.cmake
    DESTINATION "."
    RENAME TrilinosConfig.cmake
  )
  
  
  #
  # Configure the version file for Trilinos
  #
  
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/cmake/TrilinosConfigVersion.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/TrilinosConfigVersion.cmake @ONLY)

  INSTALL(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/TrilinosConfigVersion.cmake
    DESTINATION "."
  )

ENDFUNCTION()

