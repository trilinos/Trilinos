################################################################################
#
# Add target for gathering up build stats
#
################################################################################


include("${CMAKE_CURRENT_LIST_DIR}/BuildStatsSharedVars.cmake")


# Create custom 'gather-build-stats' target that will run last
#
# NOTE: This function must be called at the very end of all of the build
# targets that get created for a project!
#
function(add_target_gather_build_stats)

  if (${PROJECT_NAME}_ENABLE_BUILD_STATS)

    add_custom_command(
      OUTPUT "${BUILD_STATS_CSV_FILE}"
      COMMAND "${BUILD_STATS_SRC_DIR}/gather_build_stats.py"
      WORKING_DIRECTORY "${${PROJECT_NAME}_BINARY_DIR}" )

    add_custom_target(gather-build-stats ALL
      DEPENDS "${BUILD_STATS_CSV_FILE}")

    get_all_build_targets_including_in_subdirs("${${PROJECT_NAME}_SOURCE_DIR}"
      projectBuildTargetsList)

    if (projectBuildTargetsList)
      add_dependencies(gather-build-stats ${projectBuildTargetsList})
    endif()

  endif()

endfunction()


# Get a list all of the lib and exec build targets starting in a subdir and in
# below subdirs.
#
function(get_all_build_targets_including_in_subdirs  srcdir  targetsListVarOut)

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
