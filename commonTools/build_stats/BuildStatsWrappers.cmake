################################################################################
#
# Set up for build stats compiler wrappers and gathering up build stats
#
################################################################################

set(BUILD_STATS_SRC_DIR "${CMAKE_CURRENT_LIST_DIR}")

set(BUILD_STATS_CSV_FILE "${${PROJECT_NAME}_BINARY_DIR}/build_stats.csv")

function(build_stats_find_time)

  # let the user provide BUILD_STATS_TIME_CMD
  if (BUILD_STATS_TIME_CMD)
    message("-- BUILD_STATS_TIME_CMD=${BUILD_STATS_TIME_CMD}")
    set(GNU_TIME_EXE "${BUILD_STATS_TIME_CMD}")
  else()
    find_program(GNU_TIME_EXE "time" HINTS "/usr/bin")
    if(GNU_TIME_EXE)
      message("-- Found time at ${GNU_TIME_EXE}")
    else()
      message("-- GNU time NOT found")
      message("-- Install GNU time and/or set BUILD_STATS_TIME_CMD=/path/to/time")
      return()
    endif()
  endif()

  # this should ideally call the python script and request the fields
  # to test, add 'badflag' or some other nonsense
  SET(GNU_TIME_POSSIBLE_FIELDS "e;M;K;D;X;F;R;W;w;c;S;U;P;I;O;r;s;k;x")
  SET(GNU_TIME_SUPPORTED_FIELDS "")
  
  #  should ideally ask for the dtypes or suitable regexes to vet them
  foreach(flag ${GNU_TIME_POSSIBLE_FIELDS})
    message(DEBUG "----------------------")
    message(DEBUG "Time: Testing field ${flag}")
    # the output from time goes to stderr, the programs output to stdout
    set(GNU_TIME_RC "-1")
    set(GNU_TIME_OUTPUT "")
    execute_process(COMMAND "${GNU_TIME_EXE}"
                              "--format=%${flag}" "true"
                    # this is useless - we run a noop command
                    RESULT_VARIABLE GNU_TIME_RC
                    # capture stderr
                    ERROR_VARIABLE GNU_TIME_OUTPUT
                    )
    # if this fails, then something is broken on the system.
    # the checks after will likely fail, because they expect a predefined
    # format for stderr text
    if(GNU_TIME_RC MATCHES "0")
      #message("Time invocation error returned `${GNU_TIME_RC}` but expected `0`")
    else()
      message(DEBUG "Time invocation error returned `${GNU_TIME_RC}` but expected `0`")
      message("-- GNU_TIME_EXE=${GNU_TIME_EXE} does not work")
      message("-- Set BUILD_STATS_TIME_CMD if /usr/bin/time is invalid")
      unset(BUILD_STATS_TIME_CMD CACHE)
      return()
    endif()

    # this isn't ideal, we match the failure mode, but we should match valid output
    # That is harder to do, as we would need to determine expected format for each
    # command (e.g., dtype or a suitable regex from the python script)
    # A good todo would be to have the pythong script report headeres it will
    # try and a valid expected format regex, then we test against that here.
    if("${GNU_TIME_OUTPUT}" MATCHES "^?${flag}.*")
      message(DEBUG "Time does not support Field: ${flag}")
      # not ideal, but let's quit if we don't find the fields...
      # a better solution would be to disable the field, and only
      # error out if we can't find some required set of fields.
      message("-- GNU_TIME_EXE=${GNU_TIME_EXE} does not work")
      message("-- Set BUILD_STATS_TIME_CMD if /usr/bin/time is invalid")
      unset(BUILD_STATS_TIME_CMD CACHE)
      return()
    else()
      message(DEBUG "Time supports Field: ${flag}")
      list(APPEND GNU_TIME_SUPPORTED_FIELDS "${flag}")
    endif()
  endforeach()

  # if we get here, we should have a list of supported fields from TIME
  # We could reconcile the supported fields against anything the user specified
  # in the TRILINOS_BUILD_STATS_OUTPUT_FIELDS, which will override the default
  # fields the parser will output
  # `time` will actually silently accept bad fields, and give `?field` back
  # if we use TRILINOS_BUILD_STATS_OUTPUT_FIELDS then bad fields will simply not
  # be written to a file
  #
  # an unimplemented feature in the wrapper is `TRILINOS_BUILD_STATS_PARSE_NM`
  # which could control if NM is used. (like time, we expect it to work)
  # we could `find_program` it as well

  set(BUILD_STATS_TIME_CMD ${GNU_TIME_EXE}
      CACHE FILEPATH "The GNU time binary required by build_stats" FORCE )
endfunction()


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
    # sets BUILD_STATS_TIME_CMD
    build_stats_find_time()

    # give up if this didn't work
    if(NOT BUILD_STATS_TIME_CMD)
      message("-- ${PROJECT_NAME}_ENABLE_BUILD_STATS=ON, but valid GNU Time was not found")
      message("-- Disabling BUILD_STATS")
      unset(${PROJECT_NAME}_ENABLE_BUILD_STATS CACHE)
      return()
    endif()

    generate_build_stats_wrapper_for_op(C   WRAP CMAKE_C_COMPILER)
    generate_build_stats_wrapper_for_op(CXX WRAP CMAKE_CXX_COMPILER)
    if (${PROJECT_NAME}_ENABLE_Fortran)
      generate_build_stats_wrapper_for_op(Fortran WRAP CMAKE_Fortran_COMPILER)
    endif()

    generate_build_stats_wrapper_for_op(LD WRAP CMAKE_LD ALLOW_FIND)
    generate_build_stats_wrapper_for_op(AR WRAP CMAKE_AR ALLOW_FIND)
    generate_build_stats_wrapper_for_op(RANLIB WRAP CMAKE_RANLIB ALLOW_FIND)

  endif()

endfunction()


# Generate the build stats compiler wrapper for a given CMake variable.
#
# The intent of this function is pass in arbitrary cmake variables that
# map to commands and generate suitable wrappers.
#
# Supported functions are C, CXX, Fortran, AR, LD, and RANLIB
#
#  TODO: wrap MPIEXEC!
function(generate_build_stats_wrapper_for_op op_name)
  cmake_parse_arguments(
        PARSE_ARGV 1
        BUILD_STATS
        "ALLOW_FIND"
        "WRAP"
        ""
  )
  if (NOT BUILD_STATS_TIME_CMD)
    message(DEBUG "generate_build_stats_wrapper_for_op: BUILD_STATS_TIME_CMD is required but it is not set")
    return()
  endif()

  set(variable_to_set ${BUILD_STATS_WRAP})
  set(ALLOW_FIND ${BUILD_STATS_ALLOW_FIND})

  # 'op_name' is the short name, like C, CXX, Fortran, LD, AR, RANLIB
  # this is will give us a lowercase 
  string(TOLOWER "${op_name}" op_lc)
  set(op_wrapper
    "${${PROJECT_NAME}_BINARY_DIR}/build_stat_${op_lc}_wrapper.sh")

  # there's an issue here - if CMAKE_FOO is unset (whatever `variable_to_set` is)
  # we need a to know the command - but CMake hasn't chosen one yet...
  if( (ALLOW_FIND) AND ( ("${${variable_to_set}}" STREQUAL "")
                         OR
                         (NOT ${variable_to_set}) ))
    message("-- " "${variable_to_set} is not set, but a wrapper has been requested. Asking CMake to find ${op_lc}")
    find_program(${variable_to_set} "${op_lc}")
    print_var(${variable_to_set})
  endif()

  # Override the op with the wrapper but remember the original command
  # we take as a parameter a CMake `variable_to_set`, e.g., CMAKE_CXX_COMPILER
  # we only wrap if `variable_to_set` is defined.
  if( ("${${variable_to_set}_ORIG}" STREQUAL "") AND (${variable_to_set}) )
    # we want to set CMAKE_FOO_ORIG to CMAKE_FOO (e.g., FOO=C_COMPILER or AR)
    # `variable_to_set` is a string, so we need the value of the string evaluated
    set(${variable_to_set}_ORIG ${${variable_to_set}}
      CACHE FILEPATH "Original non-wrappeed ${op_name}" FORCE )
    set(${variable_to_set} "${op_wrapper}"
      CACHE FILEPATH "Overwritten build stats ${op_name} wrapper" FORCE )

    message("-- " "Generating build stats wrapper for ${op_name}")
    # this variable is only used to write the subsequent configure file
    # so the variable name is generic.
    set(BUILD_STATS_WRAPPER_INNER_OP "${${variable_to_set}_ORIG}")
    configure_file("${BUILD_STATS_SRC_DIR}/build_stat_lang_wrapper.sh.in"
                   "${op_wrapper}" @ONLY)

    # Use the orginal compiler for the installed <XXX>Config.cmake files
    # doubt this works w/AR/LD/RANLIB
    set(${variable_to_set}_OP_FOR_CONFIG_FILE_INSTALL_DIR
        "${${variable_to_set}_ORIG}" CACHE INTERNAL "")

  # conditional expanded so we only wrap what is defined
  else()
     # allow some verbosity when we don't set something
     message("-- " "Not wrapping ${op_name} because "
                     "${variable_to_set}=`${variable_to_set}` is not set."
                     " To enable statistics set ${variable_to_set}.")
  endif()
endfunction()
# NOTE: The above implementation will make sure the compiler wrapper will get
# updated if the *.sh.in template file changes and just reconfiguring.
# Actaully, you should be able to fix the wrapper and just type 'make' and it
# should reconfigure and update automatically.


# Get the var BASE_BUILD_DIR_FOR_PYTHON
#
macro(get_base_build_dir_for_python)
  set(get_cwd_for_python ${BUILD_STATS_SRC_DIR}/get_cwd_for_python.py)
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} ${get_cwd_for_python}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    OUTPUT_VARIABLE BASE_BUILD_DIR_FOR_PYTHON
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  print_var(BASE_BUILD_DIR_FOR_PYTHON)
endmacro()
# NOTE: We need this function to get the value of os.getcwd() from Python so
# that it matches the value returned inside of magic_wapper.py.  The issue is
# that some platforms, CMake determines a different absolute base build dir
# for systems with mounted filesystems.  The only systems I know this happens
# on are some systems at SNL with the mounted home directories.  By using the
# same Python code, we ensure that we get the same base directory, which is
# needed when computing relative paths.


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

  # disable this for now...
  return()

  if (${PROJECT_NAME}_ENABLE_BUILD_STATS)
    install_build_stats_wrapper_for_lang(C)
    install_build_stats_wrapper_for_lang(CXX)
    if (${PROJECT_NAME}_ENABLE_Fortran)
      install_build_stats_wrapper_for_lang(Fortran)
    endif()

    install_build_stats_wrapper_for_lang(AR)
    install_build_stats_wrapper_for_lang(LD)
    install_build_stats_wrapper_for_lang(RANLIB)
  endif()

endfunction()


# Install the build stats compiler wrapper for a single language.
#
function(install_build_stats_wrapper_for_lang  lang)
  string(TOLOWER "${op_name}" op_lc)
  set(op_wrapper
    "${${PROJECT_NAME}_BINARY_DIR}/build_stat_${op_lc}_wrapper.sh")

  install(PROGRAMS "${op_wrapper}"
    DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}")
endfunction()


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
