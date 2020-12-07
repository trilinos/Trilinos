################################################################################
#
# Set up for build stats compiler wrappers and gathering up build stats
#
################################################################################

set(BUILD_STATS_SRC_DIR "${CMAKE_CURRENT_LIST_DIR}")

set(BUILD_STATS_CSV_FILE "${${PROJECT_NAME}_BINARY_DIR}/build_stats.csv")

# Generate the build stats compiler wrappers if asked to do so.
#
function(generate_build_stats_wrappers)

  # Set default for cache var ${PROJECT_NAME}_ENABLE_BUILD_STATS
  if (NOT "$ENV{${PROJECT_NAME}_ENABLE_BUILD_STATS}" STREQUAL "")
    # Use the default set in the env (overrides any local default set)
    set(${PROJECT_NAME}_ENABLE_BUILD_STATS_DEFAULT
      "$ENV{${PROJECT_NAME}_ENABLE_BUILD_STATS}")
  elseif(NOT "${${PROJECT_NAME}_ENABLE_BUILD_STATS_DEFAULT}" STREQUAL "")
    # ${PROJECT_NAME}_ENABLE_BUILD_STATS_DEFAULT was already set, so use it as
    # the default.
  else()
    # No default was set, so make it OFF by default
    set(${PROJECT_NAME}_ENABLE_BUILD_STATS_DEFAULT OFF)
  endif()

  # Set cache var ${PROJECT_NAME}_ENABLE_BUILD_STATS
  advanced_set(${PROJECT_NAME}_ENABLE_BUILD_STATS
    ${${PROJECT_NAME}_ENABLE_BUILD_STATS_DEFAULT} CACHE BOOL
    "If set to 'ON', then compiler wrappers will be created and used to gather build stats."
    )

  # Generate the build-stats compiler wrappers
  get_base_build_dir_for_python()
  if (${PROJECT_NAME}_ENABLE_BUILD_STATS)
    generate_build_stats_wrapper_for_lang(C)
    generate_build_stats_wrapper_for_lang(CXX)
    if (${PROJECT_NAME}_ENABLE_Fortran)
      generate_build_stats_wrapper_for_lang(Fortran)
    endif()

    set(gather_build_status "${${PROJECT_NAME}_BINARY_DIR}/gather_build_stats.sh")
    configure_file("${BUILD_STATS_SRC_DIR}/gather_build_stats.sh"
      "${gather_build_status}" COPYONLY)
  endif()

endfunction()


# Get the var BASE_BUILD_DIR_FOR_PYTHON
#
macro(get_base_build_dir_for_python)
  set(get_cwd_for_python ${BUILD_STATS_SRC_DIR}/get_cwd_for_python.py)
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} ${get_cwd_for_python}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    OUTPUT_VARIABLE BASE_BUILD_DIR_FOR_PYTHON
    OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro()
# NOTE: We need this function to get the value of os.getcwd() from Python so
# that it matches the value returned inside of magic_wapper.py.  The issue is
# that some platforms, CMake determines a different absolute base build dir
# for systems with mounted filesystems.  The only systems I know this happens
# on are some systems at SNL with the mounted home directories.  By using the
# same Python code, we ensure that we get the same base directory, which is
# needed when computing relative paths.


# Generate the build stats compiler wrapper for a single language.
#
function(generate_build_stats_wrapper_for_lang  lang)

  string(TOLOWER "${lang}" lang_lc)
  set(compiler_wrapper
    "${${PROJECT_NAME}_BINARY_DIR}/build_stat_${lang_lc}_wrapper.sh")

  # Override the compiler with the wrapper but remember the original compiler
  if ("${CMAKE_${lang}_COMPILER_ORIG}" STREQUAL "")
    set(CMAKE_${lang}_COMPILER_ORIG "${CMAKE_${lang}_COMPILER}"
      CACHE FILEPATH "Original non-wrappeed ${lang} compiler" FORCE )
    set(CMAKE_${lang}_COMPILER "${compiler_wrapper}"
      CACHE FILEPATH "Overwritten build stats ${lang} compiler wrapper" FORCE )
  endif()

  message("-- " "Generate build stats compiler wrapper for ${lang}")
  set(BUILD_STAT_COMPILER_WRAPPER_INNER_COMPILER "${CMAKE_${lang}_COMPILER_ORIG}")
  configure_file("${BUILD_STATS_SRC_DIR}/build_stat_lang_wrapper.sh.in"
    "${compiler_wrapper}" @ONLY)

  # Use the orginal compiler for the installed <XXX>Config.cmake files
  set(CMAKE_${lang}_COMPILER_FOR_CONFIG_FILE_INSTALL_DIR
    "${CMAKE_${lang}_COMPILER_ORIG}" CACHE INTERNAL "")
  #print_var(CMAKE_${lang}_COMPILER_FOR_CONFIG_FILE_INSTALL_DIR)

endfunction()

# NOTE: The above implementation will make sure the compiler wrapper will get
# updated if the *.sh.in template file changes and just reconfiguring.
# Actaully, you should be able to fix the wrapper and just type 'make' and it
# should reconfigure and update automatically.


# Remove the build stats file on configure if asked to do so.
#
function(remove_build_stats_file_on_configure)

  advanced_set(${PROJECT_NAME}_REMOVE_BUILD_STATS_ON_CONFIGURE OFF
    ${${PROJECT_NAME}_REMOVE_BUILD_STATS_ON_CONFIGURE_DEFAULT} CACHE BOOL
    "If set to 'ON', then the build_stats.csv file will be removed on each configure."
    )

  if (
      (${PROJECT_NAME}_REMOVE_BUILD_STATS_ON_CONFIGURE)
      AND
      (EXISTS "${BUILD_STATS_CSV_FILE}")
    )
    MESSAGE("-- " "Removing existing file '${BUILD_STATS_CSV_FILE}'")
    file(REMOVE "${BUILD_STATS_CSV_FILE}")
  endif()

endfunction()


# Remove the <target>.timing files on a fresh configure if asked to do so.
#
function(remove_build_stats_timing_files_on_fresh_configure)

  advanced_set(${PROJECT_NAME}_REMOVE_BUILD_STATS_TIMING_FILES_ON_FRESH_CONFIGURE OFF
    CACHE BOOL
    "If set to 'ON', then all <target>.timing files will be removed on a freash configure."
    )

  if (
    ${PROJECT_NAME}_REMOVE_BUILD_STATS_TIMING_FILES_ON_FRESH_CONFIGURE
    AND
    (NOT "${${PROJECT_NAME}_BUILD_STATS_INIT_CONFIG_WAS_DONE}")
    )

    message("-- " "Removing all <target>.timing files on fresh configure")

    execute_process(
      COMMAND "${BUILD_STATS_SRC_DIR}/remove_all_target_timing_files.sh"
        ${PROJECT_BINARY_DIR} )

    set(${PROJECT_NAME}_BUILD_STATS_INIT_CONFIG_WAS_DONE ON CACHE INTERNAL "")

  endif()

endfunction()


# Set up install targets for the build stats scripts
#
# This can only be called after these install dirs have been set!
#
function(install_build_stats_scripts)

  if (${PROJECT_NAME}_ENABLE_BUILD_STATS)

    set(gather_build_status "${${PROJECT_NAME}_BINARY_DIR}/gather_build_stats.sh")
    install(PROGRAMS "${gather_build_status}" 
      DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}")

    install_build_stats_wrapper_for_lang(C)
    install_build_stats_wrapper_for_lang(CXX)
    if (${PROJECT_NAME}_ENABLE_Fortran)
      install_build_stats_wrapper_for_lang(Fortran)
    endif()

  endif()

endfunction()


# Install the build stats compiler wrapper for a single language.
#
function(install_build_stats_wrapper_for_lang  lang)
  string(TOLOWER "${lang}" lang_lc)
  set(compiler_wrapper
    "${${PROJECT_NAME}_BINARY_DIR}/build_stat_${lang_lc}_wrapper.sh")
  install(PROGRAMS "${compiler_wrapper}"
    DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}")
endfunction()


# Create custom 'gather-build-stats' target that will run last
#
# NOTE: This function must be called at the very end of all of the build
# targets that get created for a project!
#
function(add_target_gather_build_stats)

  if (${PROJECT_NAME}_ENABLE_BUILD_STATS)

    set(buildStatsCsvFile "${${PROJECT_NAME}_BINARY_DIR}/build_stats.csv")

    add_custom_command(
      OUTPUT "${buildStatsCsvFile}"
      COMMAND "${${PROJECT_NAME}_BINARY_DIR}/gather_build_stats.sh"
      WORKING_DIRECTORY "${${PROJECT_NAME}_BINARY_DIR}" )

    add_custom_target(gather-build-stats ALL
      DEPENDS "${buildStatsCsvFile}")

    get_all_build_targets_including_in_subdirs("${${PROJECT_NAME}_SOURCE_DIR}"
      projectBuildTargetsList)

    if (projectBuildTargetsList)
      add_dependencies(gather-build-stats ${projectBuildTargetsList})
    endif()

  endif()

endfunction()


# Get a list all of the lib and exec build targets starting in a a subdir and
# in below subdirs.
#
function(get_all_build_targets_including_in_subdirs srcdir  targetsListVarOut)

  set(targetsList "")

  # Recurse into subdirectories.
  get_property(dirs DIRECTORY ${srcdir} PROPERTY SUBDIRECTORIES)
  foreach(d IN LISTS dirs)
    get_all_build_targets_including_in_subdirs(${d} targetsSubdirList)
    list(APPEND targetsList ${targetsSubdirList})
  endforeach()

  # Get the targets from this directory.
  get_property(allTargetsThisDir DIRECTORY ${srcdir} PROPERTY BUILDSYSTEM_TARGETS)
  filter_only_build_targets(allTargetsThisDir buildTargetsThisDir)
  list(APPEND targetsList ${buildTargetsThisDir})

  # Return
  set(${targetsListVarOut} ${targetsList} PARENT_SCOPE)

endfunction()


function(filter_only_build_targets targetListInVar targetListOutVar)

  #print_var(targetListInVar)
  #print_var(${targetListInVar})

  set(targetListOut "")

  foreach (target IN LISTS ${targetListInVar})
    #print_var(target)
    get_property(targetType TARGET ${target} PROPERTY TYPE)
    #print_var(targetType)
    if (
        targetType STREQUAL "STATIC_LIBRARY" OR
        targetType STREQUAL "SHARED_LIBRARY" OR
        targetType STREQUAL "EXECUTABLE"
      )
      #message("-- " "${target} is a regular build target!")
      list(APPEND targetListOut ${target})
    else()
      #message("-- " "${target} is **NOT** a regular build target!")
    endif()
  endforeach()

  set(${targetListOutVar} ${targetListOut} PARENT_SCOPE)

endfunction()
