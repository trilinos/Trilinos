# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# TriBITS package_arch includes
include(TribitsConfigureTiming)
include(TribitsGetPackageSublists)

# TriBITS utils includes
include(FindListElement)


# @MACRO: tribits_exclude_files()
#
# Exclude package files/dirs from the source distribution by appending
# ``CPACK_SOURCE_IGNORE_FILES``.
#
# Usage::
#
#  tribits_exclude_files(<file0> <file1> ...)
#
# This is called in the top-level parent package's
# `<packageDir>/CMakeLists.txt`_ file and each file or directory name
# ``<filei>`` is actually interpreted by CMake/CPack as a regex that is
# prefixed by the project's and package's source directory names so as to not
# exclude files and directories of the same name and path from other packages.
# If ``<filei>`` is an absolute path it is not prefixed but is appended to
# ``CPACK_SOURCE_IGNORE_FILES`` unmodified.
#
# In general, do **NOT** put in excludes for files and directories that are
# not under this package's source tree.  If the given package is not enabled,
# then this command will never be called! For example, don't put in excludes
# for PackageB's files in PackageA's ``CMakeLists.txt`` file because if
# PackageB is enabled but PackageA is not, the excludes for PackageB will
# never get added to ``CPACK_SOURCE_IGNORE_FILES``.
#
# Also, be careful to note that the ``<filei>`` arguments are actually regexes
# and one must be very careful to understand how CPack will use these regexes
# to match files that get excluded from the tarball.  For more details, see
# `Creating Source Distributions`_.
#
macro(tribits_exclude_files)

  if (NOT "${${PACKAGE_NAME}_PARENT_PACKAGE}" STREQUAL "")
    message(FATAL_ERROR
      "ERROR: tribits_exclude_files() was called in a subpackage CmakeLists.txt file!"
      "  Instead, move this call to the file"
      " ${${${PACKAGE_NAME}_PARENT_PACKAGE}_SOURCE_DIR}/CMakeLists.txt"
      " and adjust the paths accordingly!" )
  endif()

  set(FILES_TO_EXCLUDE ${ARGN})

  # Need to add "/<project source dir>/<package dir>/" to each file to prevent
  # someone from trying to exclude a file like "readme" and having it
  # inadvertently exclude a file matching that name in another package.
  set(MODIFIED_FILES_TO_EXCLUDE "")

  set(${PROJECT_NAME}_SOURCE_PATH ${${PROJECT_NAME}_SOURCE_DIR})

  foreach(FILE ${FILES_TO_EXCLUDE})
    #Ensure that if the full path was specified for the file that we don't add
    #"/<project source dir>/<package dir>/" again.
    set(MATCH_STRING "${${PACKAGE_NAME}_SOURCE_DIR}")
    string(REGEX MATCH ${MATCH_STRING} MATCHED ${FILE} )
    if(NOT MATCHED)
      list(APPEND MODIFIED_FILES_TO_EXCLUDE
        "${${PACKAGE_NAME}_SOURCE_DIR}/${FILE}")
    else()
      list(APPEND MODIFIED_FILES_TO_EXCLUDE ${FILE})
    endif()
  endforeach()

#Leaving in for debugging purposes
#  message("List of files being excluded for package ${PACKAGE_NAME}")
#  foreach(NEW_FILE ${MODIFIED_FILES_TO_EXCLUDE})
#    message(${NEW_FILE})
#  endforeach()

  list(APPEND CPACK_SOURCE_IGNORE_FILES ${MODIFIED_FILES_TO_EXCLUDE})
  if (NOT ${PROJECT_NAME}_BINARY_DIR STREQUAL ${PACKAGE_NAME}_BINARY_DIR)
    set(CPACK_SOURCE_IGNORE_FILES ${CPACK_SOURCE_IGNORE_FILES} PARENT_SCOPE)
  endif()

endmacro()


# Set up for packaging and distribution
#
macro(tribits_setup_packaging_and_distribution)

  tribits_config_code_start_timer(CPACK_SETUP_TIME_START_SECONDS)

  # K.1) Run callback function for the base project.

  tribits_project_define_packaging_runner()
  # The above must define the basic project settings for CPACK that are
  # specific to the project and should not be provided by the user.

  # K.2) Removing any packages or packages not enabled from the tarball

  if (${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION)
    set(tribitsPackageList ${${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES})
  else()
    set(tribitsPackageList ${${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES})
  endif()

  tribits_get_sublist_nonenabled(tribitsPackageList  nonEnabledTribitsPackage  "")

  foreach(TRIBITS_PACKAGE ${nonEnabledTribitsPackage})

    # Determine if this is a package to not ignore
    find_list_element(TRIBITS_CPACK_PACKAGES_TO_NOT_IGNORE
       ${TRIBITS_PACKAGE}  TRIBITS_PACKAGE_DONT_IGNORE)

    if (NOT TRIBITS_PACKAGE_DONT_IGNORE)

      # Checking if we have a relative path to the package's files. Since the
      # exclude is a regular expression any "../" will be interpreted as <any
      # char><any char>/ which would never match the package's actual
      # directory. There isn't a direct way in cmake to convert a relative
      # path into an absolute path with string operations so as a way of
      # making sure that we get the correct path of the package we use a
      # find_path for the CMakeLists.txt file for the package. Since the
      # package has to have this file to work correctly it should be
      # guaranteed to be there.
      string(REGEX MATCH "[.][.]/" RELATIVE_PATH_CHARS_MATCH
        ${${TRIBITS_PACKAGE}_REL_SOURCE_DIR})
      if ("${RELATIVE_PATH_CHARS_MATCH}" STREQUAL "")
        list(PREPEND CPACK_SOURCE_IGNORE_FILES
          "${PROJECT_SOURCE_DIR}/${${TRIBITS_PACKAGE}_REL_SOURCE_DIR}/")
      else()
        find_path(ABSOLUTE_PATH  CMakeLists.txt  PATHS
          "${PROJECT_SOURCE_DIR}/${${TRIBITS_PACKAGE}_REL_SOURCE_DIR}"
          NO_DEFAULT_PATH)
        if ("${ABSOLUTE_PATH}" STREQUAL "ABSOLUTE_PATH-NOTFOUND")
          message(AUTHOR_WARNING "Relative path found for disabled package"
            " ${TRIBITS_PACKAGE} but package was missing a CMakeLists.txt file."
            " This disabled package will likely not be excluded from a source release")
        endif()
        list(PREPEND CPACK_SOURCE_IGNORE_FILES "${ABSOLUTE_PATH}")
      endif()
    endif()

  endforeach()

  # Add excludes for VC files/dirs
  list(APPEND CPACK_SOURCE_IGNORE_FILES
    /[.]git/
    [.]gitignore$
    )

  # K.3) Set up install component dependencies

  tribits_get_sublist_enabled(
    ${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES
    enabledInternalToplevelPackages  "")

  foreach(pkgName ${enabledInternalToplevelPackages})
    if(NOT "${${pkgName}_LIB_ENABLED_DEPENDENCIES}" STREQUAL "")
      string(TOUPPER ${pkgName} upperPkgName)
      set(CPACK_COMPONENT_${upperPkgName}_DEPENDS ${${pkgName}_LIB_ENABLED_DEPENDENCIES})
      # ToDo: The above needs to be changed to the list of *internal* enabled
      # package dependencies!  (But there are no tests for this currently and
      # I am not sure who is using this.)
    endif()
  endforeach()

  # K.4) Resetting the name to avoid overwriting registry keys when installing

  if(WIN32)
    set(CPACK_PACKAGE_NAME "${CPACK_PACKAGE_NAME}-${${PROJECT_NAME}_VERSION}")
    if (TPL_ENABLE_MPI)
      set(CPACK_PACKAGE_NAME "${CPACK_PACKAGE_NAME}-mpi")
    ELSE ()
      set(CPACK_PACKAGE_NAME "${CPACK_PACKAGE_NAME}-serial")
    endif()
    set(CPACK_GENERATOR "NSIS")
    set(CPACK_NSIS_MODIFY_PATH OFF)
  endif()

  # K.5) Determine the source generator
  if ("${${PROJECT_NAME}_CPACK_SOURCE_GENERATOR_DEFAULT}" STREQUAL "")
    set(${PROJECT_NAME}_CPACK_SOURCE_GENERATOR_DEFAULT "TGZ")
  endif()
  set(${PROJECT_NAME}_CPACK_SOURCE_GENERATOR
    ${${PROJECT_NAME}_CPACK_SOURCE_GENERATOR_DEFAULT}
    CACHE STRING
    "The types of source generators to use for CPACK_SOURCE_GENERATOR.")
  set(CPACK_SOURCE_GENERATOR ${${PROJECT_NAME}_CPACK_SOURCE_GENERATOR})

  # K.6) Loop through the Repositories and run their callback functions.
  foreach(REPO ${${PROJECT_NAME}_ALL_REPOSITORIES})
    tribits_get_repo_name_dir(${REPO}  REPO_NAME  REPO_DIR)
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("Processing packaging call-backs for ${REPO_NAME}")
    endif()
    tribits_repository_define_packaging_runner(${REPO_NAME})
  endforeach()

  # K.7) Include <Project>RepoVersion.txt if generated
  set(PROJECT_REPO_VERSION_FILE
     "${CMAKE_CURRENT_BINARY_DIR}/${${PROJECT_NAME}_REPO_VERSION_FILE_NAME}")
  if (EXISTS "${PROJECT_REPO_VERSION_FILE}")
    foreach(SOURCE_GEN ${CPACK_SOURCE_GENERATOR})
      set(CPACK_INSTALL_COMMANDS ${CPACK_INSTALL_COMMANDS}
        "${CMAKE_COMMAND} -E copy '${PROJECT_REPO_VERSION_FILE}' '${CMAKE_CURRENT_BINARY_DIR}/_CPack_Packages/Linux-Source/${SOURCE_GEN}/${CPACK_PACKAGE_NAME}-${${PROJECT_NAME}_VERSION}-Source/${${PROJECT_NAME}_REPO_VERSION_FILE_NAME}'")
    endforeach()
  endif()

  # Print the set of excluded files
  if(${PROJECT_NAME}_VERBOSE_CONFIGURE OR
    ${PROJECT_NAME}_DUMP_CPACK_SOURCE_IGNORE_FILES
    )
    message("Exclude files when building source packages:")
    foreach(item IN LISTS CPACK_SOURCE_IGNORE_FILES)
      message(${item})
    endforeach()
  endif()

  # K.8) Finally process with CPack
  include(CPack)

  tribits_config_code_stop_timer(CPACK_SETUP_TIME_START_SECONDS
    "Total time to set up for CPack packaging")

endmacro()


macro(tribits_project_define_packaging_runner)
  set(CALLBACK_DEFINE_PACKAGING_FILE
    "${PROJECT_SOURCE_DIR}/cmake/CallbackDefineProjectPackaging.cmake")
  #print_var(CALLBACK_DEFINE_PACKAGING_FILE)
  if (EXISTS ${CALLBACK_DEFINE_PACKAGING_FILE})
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("Processing call-back file and macros in"
        " '${CALLBACK_DEFINE_PACKAGING_FILE}'")
    endif()
    # Define the callback macros as empty in case it is not defined
    # in this file.
    create_empty_tribits_project_define_packaging()
    # Include the file which will define the callback macros
    tribits_trace_file_processing(PROJECT  INCLUDE
      "${CALLBACK_DEFINE_PACKAGING_FILE}")
    include(${CALLBACK_DEFINE_PACKAGING_FILE})
    # Call the callback macros to inject project-specific behavir
    tribits_project_define_packaging()
    # Set back the callback macros to empty to ensure that no-one calls them
    create_empty_tribits_project_define_packaging()
  endif()
endmacro()


macro(create_empty_tribits_repository_setup_extra_options)
  macro(tribits_repository_setup_extra_options)
  endmacro()
endmacro()


macro(tribits_repository_setup_extra_options_runner  REPO_NAME)
  set(CALLBACK_SETUP_EXTRA_OPTIONS_FILE
    "${${REPO_NAME}_SOURCE_DIR}/cmake/CallbackSetupExtraOptions.cmake")
  #print_var(CALLBACK_SETUP_EXTRA_OPTIONS_FILE)
  if (EXISTS ${CALLBACK_SETUP_EXTRA_OPTIONS_FILE})
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("Processing call-back file and macros in"
        " '${CALLBACK_SETUP_EXTRA_OPTIONS_FILE}'")
    endif()
    # Define the callback macros as empty in case it is not defined
    # in this file.
    create_empty_tribits_repository_setup_extra_options()
    # Include the file which will define the callback macros
    set(REPOSITORY_NAME ${REPO_NAME})
    tribits_trace_file_processing(REPOSITORY  INCLUDE
      "${CALLBACK_SETUP_EXTRA_OPTIONS_FILE}")
    include(${CALLBACK_SETUP_EXTRA_OPTIONS_FILE})
    # Call the callback macros to inject repository-specific behavir
    tribits_repository_setup_extra_options()
    # Set back the callback macros to empty to ensure that nonone calls them
    create_empty_tribits_repository_setup_extra_options()
  endif()
endmacro()


macro(create_empty_tribits_repository_define_packaging)
  macro(tribits_repository_define_packaging)
  endmacro()
endmacro()


macro(tribits_repository_define_packaging_runner  REPO_NAME)
  set(CALLBACK_DEFINE_PACKAGING_FILE
    "${${REPO_NAME}_SOURCE_DIR}/cmake/CallbackDefineRepositoryPackaging.cmake")
  #print_var(CALLBACK_DEFINE_PACKAGING_FILE)
  if (EXISTS ${CALLBACK_DEFINE_PACKAGING_FILE})
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("Processing call-back file and macros in"
        " '${CALLBACK_DEFINE_PACKAGING_FILE}'")
    endif()
    # Define the callback macros as empty in case it is not defined
    # in this file.
    create_empty_tribits_repository_define_packaging()
    # Include the file which will define the callback macros
    tribits_trace_file_processing(REPOSITORY  INCLUDE
      "${CALLBACK_DEFINE_PACKAGING_FILE}")
    include(${CALLBACK_DEFINE_PACKAGING_FILE})
    # Call the callback macros to inject repository-specific behavir
    tribits_repository_define_packaging()
    # Set back the callback macros to empty to ensure that nonone calls them
    create_empty_tribits_repository_define_packaging()
  endif()
endmacro()


macro(create_empty_tribits_project_define_packaging)
  macro(tribits_project_define_packaging)
  endmacro()
endmacro()
