# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


#
# Find an installed version of ${PROJECT_NAME} for installation testing
# (the check that we are in installation mode is inside the macro)
#

function(tribits_find_project_install)

  if(${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
    if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("Searching for ${PROJECT_NAME} installation at ${${PROJECT_NAME}_INSTALLATION_DIR}/include")
    endif()
    find_package(${PROJECT_NAME} REQUIRED HINTS ${${PROJECT_NAME}_INSTALLATION_DIR})

    if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("Found ${PROJECT_NAME} installation version ${${PROJECT_NAME}_VERSION} at ${${PROJECT_NAME}_DIR}")
    endif()

    #renaming some of the variables so that they do not clash with the variables of the same name.
    set(${PROJECT_NAME}_INSTALLATION_VERSION           ${${PROJECT_NAME}_VERSION}           PARENT_SCOPE)
    set(${PROJECT_NAME}_INSTALLATION_INCLUDE_DIRS      ${${PROJECT_NAME}_INCLUDE_DIRS}      PARENT_SCOPE)
    set(${PROJECT_NAME}_INSTALLATION_LIBRARY_DIRS      ${${PROJECT_NAME}_LIBRARY_DIRS}      PARENT_SCOPE)
    set(${PROJECT_NAME}_INSTALLATION_LIBRARIES         ${${PROJECT_NAME}_LIBRARIES}         PARENT_SCOPE)
    set(${PROJECT_NAME}_INSTALLATION_PACKAGE_LIST      ${${PROJECT_NAME}_PACKAGE_LIST}      PARENT_SCOPE)
    set(${PROJECT_NAME}_INSTALLATION_BUILD_SHARED_LIBS ${${PROJECT_NAME}_BUILD_SHARED_LIBS} PARENT_SCOPE)
    set(${PROJECT_NAME}_INSTALLATION_TPL_INCLUDE_DIRS  ${${PROJECT_NAME}_TPL_INCLUDE_DIRS}  PARENT_SCOPE)
    set(${PROJECT_NAME}_INSTALLATION_TPL_LIBRARY_DIRS  ${${PROJECT_NAME}_TPL_LIBRARY_DIRS}  PARENT_SCOPE)
    set(${PROJECT_NAME}_INSTALLATION_TPL_LIBRARIES     ${${PROJECT_NAME}_TPL_LIBRARIES}     PARENT_SCOPE)

    foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_PACKAGE_LIST})
      set(${TRIBITS_PACKAGE}_INSTALLATION_INCLUDE_DIRS     ${${TRIBITS_PACKAGE}_INCLUDE_DIRS}     PARENT_SCOPE)
      set(${TRIBITS_PACKAGE}_INSTALLATION_LIBRARY_DIRS     ${${TRIBITS_PACKAGE}_LIBRARY_DIRS}     PARENT_SCOPE)
      set(${TRIBITS_PACKAGE}_INSTALLATION_LIBRARIES        ${${TRIBITS_PACKAGE}_LIBRARIES}        PARENT_SCOPE)
      set(${TRIBITS_PACKAGE}_INSTALLATION_TPL_INCLUDE_DIRS ${${TRIBITS_PACKAGE}_TPL_INCLUDE_DIRS} PARENT_SCOPE)
      set(${TRIBITS_PACKAGE}_INSTALLATION_TPL_LIBRARY_DIRS ${${TRIBITS_PACKAGE}_TPL_LIBRARY_DIRS} PARENT_SCOPE)
      set(${TRIBITS_PACKAGE}_INSTALLATION_TPL_LIBRARIES    ${${TRIBITS_PACKAGE}_TPL_LIBRARIES}    PARENT_SCOPE)
    endforeach()

  endif()
endfunction()
