################################################################################
#
# Set up for build stats compiler wrappers and gathering up build stats
#
################################################################################


include("${CMAKE_CURRENT_LIST_DIR}/BuildStatsSharedVars.cmake")


# Generate the build stats compiler wrappers if asked to do so.
#
function(generate_build_stats_wrappers)

  set_project_enable_build_stats_var()

  if (${PROJECT_NAME}_ENABLE_BUILD_STATS)
    build_stats_find_and_check_time()  # Sets cache var BUILD_STATS_TIME_CMD
    if(NOT BUILD_STATS_TIME_CMD)
      message("-- ${PROJECT_NAME}_ENABLE_BUILD_STATS=ON, but valid GNU Time was not found")
      message("-- NOTE: Force setting ${PROJECT_NAME}_ENABLE_BUILD_STATS=OFF!")
      set(${PROJECT_NAME}_ENABLE_BUILD_STATS OFF CACHE BOOL
        "Forced to 'OFF' since valid 'time' command not found" FORCE)
      return()
    endif()

    get_base_build_dir_for_python()

    generate_build_stats_wrapper_for_op(C   WRAP CMAKE_C_COMPILER)
    generate_build_stats_wrapper_for_op(CXX WRAP CMAKE_CXX_COMPILER)
    if (${PROJECT_NAME}_ENABLE_Fortran)
      generate_build_stats_wrapper_for_op(Fortran WRAP CMAKE_Fortran_COMPILER)
    endif()

    generate_build_stats_wrapper_for_op(LD WRAP CMAKE_LD ALLOW_FIND)
    generate_build_stats_wrapper_for_op(AR WRAP CMAKE_AR ALLOW_FIND)
    generate_build_stats_wrapper_for_op(RANLIB WRAP CMAKE_RANLIB ALLOW_FIND)
    # NOTE: LD, AR, and RANDLIB can be used even in builds where
    # BUILD_SHARED_LIBS=ON because individual add_librariy() commands can
    # request static libraries be built.

    set(BUILD_STATS_COMPLETED_FIRST_CONFIG TRUE CACHE INTERNAL "")
  endif()

endfunction()


# Macro that sets the cache var ${PROJECT_NAME}_ENABLE_BUILD_STATS
#
macro(set_project_enable_build_stats_var)

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

  advanced_set(${PROJECT_NAME}_ENABLE_BUILD_STATS
    ${${PROJECT_NAME}_ENABLE_BUILD_STATS_DEFAULT} CACHE BOOL
    "If set to 'ON', then compiler wrappers will be created and used to gather build stats."
    )

endmacro()


# Find the GNU 'time' command that is used by magic_wrapper.py to extract the
# info out of the command that it runs.
#
# If this finds the GNU 'time' command and it behaves correctly, then it sets
# the cache var BUILD_STATS_TIME_CMD on output.  If BUILD_STATS_TIME_CMD is
# already set by the user in the cache and it is found to not behave
# correctly, then BUILD_STATS_TIME_CMD will be removed from the cache.
# 
function(build_stats_find_and_check_time)

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

  # This should ideally call the python script and request the fields to test,
  # add 'badflag' or some other nonsense.
  SET(GNU_TIME_POSSIBLE_FIELDS "e;M;K;D;X;F;R;W;w;c;S;U;P;I;O;r;s;k;x")
  SET(GNU_TIME_SUPPORTED_FIELDS "")
  
  #  Should ideally ask for the dtypes or suitable regexes to vet them
  foreach(flag ${GNU_TIME_POSSIBLE_FIELDS})
    message(DEBUG "----------------------")
    message(DEBUG "Time: Testing field ${flag}")
    # The output from time goes to stderr, the programs output to stdout
    execute_process(COMMAND "${GNU_TIME_EXE}"
                              "--format=%${flag}" "true"
                    # this is useless - we run a noop command
                    RESULT_VARIABLE GNU_TIME_RC
                    # capture stderr
                    ERROR_VARIABLE GNU_TIME_OUTPUT
                    )
    # If this fails, then something is broken on the system.  The checks after
    # will likely fail, because they expect a predefined format for stderr
    # text.
    if(NOT GNU_TIME_RC EQUAL 0)
      message(DEBUG "Time invocation error returned `${GNU_TIME_RC}` but expected `0`")
      message("-- GNU_TIME_EXE=${GNU_TIME_EXE} does not work")
      message("-- Unset BUILD_STATS_TIME_CMD since '${GNU_TIME_EXE}' is invalid!")
      unset(BUILD_STATS_TIME_CMD CACHE)
      return()
    endif()

    # For now, just assert that all expected fields are supported (see
    # discussion after function of other possible options).
    if("${GNU_TIME_OUTPUT}" MATCHES "^?${flag}.*")
      message("-- Time does not support Field: ${flag}")
      message("-- GNU_TIME_EXE=${GNU_TIME_EXE} does not work")
      message("-- Unset BUILD_STATS_TIME_CMD since '${GNU_TIME_EXE}' is invalid!")
      unset(BUILD_STATS_TIME_CMD CACHE)
      return()
    else()
      message(DEBUG "-- Time supports Field: ${flag}")
      list(APPEND GNU_TIME_SUPPORTED_FIELDS "${flag}")
    endif()
  endforeach()

  # If we get here, we should have a list of supported fields from TIME.  
  set(BUILD_STATS_TIME_CMD ${GNU_TIME_EXE}
      CACHE FILEPATH "The GNU time binary required by build_stats" FORCE )
endfunction()
#
# NOTE: Above, the GNU_TIME_SUPPORTED_FIELDS list var is currently not being
# used for anything but in the future, it could be exported to the env as
# TRILINOS_BUILD_STATS_OUTPUT_FIELDS for the magic_wapper.py to use to pass in
# to the 'time' command for fields that are known to be supported.  This would
# override the default fields specified there.  Note that `time` will actually
# silently accept bad fields, and give `?field` back.  If we were to set
# TRILINOS_BUILD_STATS_OUTPUT_FIELDS the GNU_TIME_SUPPORTED_FIELDS then bad
# fields will simply not be written to a file.
#
# One unimplemented feature in the wrapper is
# `TRILINOS_BUILD_STATS_PARSE_NM` which we could control if NM is used. Like
# 'time', we expect it to work and we could `find_program()` it as well.


# Get the non-cache var BASE_BUILD_DIR_FOR_PYTHON
#
# This var gets picked up in the configure of build_stat_wrapper.sh.in.
#
macro(get_base_build_dir_for_python)
  set(get_cwd_for_python ${BUILD_STATS_SRC_DIR}/get_cwd_for_python.py)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} ${get_cwd_for_python}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    OUTPUT_VARIABLE BASE_BUILD_DIR_FOR_PYTHON
    OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro()
#
# NOTE: We need this function to get the value of os.getcwd() from Python so
# that it matches the value returned inside of magic_wapper.py.  The issue is
# that some platforms, CMake determines a different absolute base build dir
# for systems with mounted filesystems.  The only systems I know this happens
# on are some systems at SNL with the mounted home directories.  By using the
# same Python code, we ensure that we get the same base directory, which is
# needed when computing relative paths.


# Generate the build stats compiler wrapper for a given CMake variable.
#
# Usage:
#
#   generate_build_stats_wrapper_for_op(<op_name> WRAP <cmake_var> [ALLOW_FIND])
#
# The intent of this function is pass in arbitrary cmake variables <cmake_var>
# that map to commands and generate suitable wrappers.
#
# <op_name> is the short name, like C, CXX, Fortran, LD, AR, RANLIB.
#
function(generate_build_stats_wrapper_for_op  op_name)
  cmake_parse_arguments(
    PARSE_ARGV 1
    BUILD_STATS    # prefix
    "ALLOW_FIND"   # options
    "WRAP"         # one_value_keywords
    ""             # multi_value_keywords
    )
  set(variable_to_set "${BUILD_STATS_WRAP}")

  string(TOLOWER "${op_name}" op_lc)
  set(op_wrapper "${${PROJECT_NAME}_BINARY_DIR}/build_stat_${op_lc}_wrapper.sh")

  generate_build_stats_wrapper_for_op_find_op_lc() # Sets ${variable_to_set}

  # Override the op with the wrapper but remember the original command
  if (NOT BUILD_STATS_COMPLETED_FIRST_CONFIG)
    if (${variable_to_set}) # True if was set on input or was found above
      set(${variable_to_set}_ORIG ${${variable_to_set}}
        CACHE FILEPATH "Original non-wrapped ${op_name}" FORCE )
      set(${variable_to_set} "${op_wrapper}"
        CACHE FILEPATH "Overwritten build stats ${op_name} wrapper" FORCE )

      message("-- " "Generating build stats wrapper for ${op_name}")
      set(BUILD_STATS_WRAPPER_INNER_OP "${${variable_to_set}_ORIG}")
      configure_file("${BUILD_STATS_SRC_DIR}/build_stat_wrapper.sh.in"
        "${op_wrapper}" @ONLY)

      set(${variable_to_set}_OP_FOR_CONFIG_FILE_INSTALL_DIR
        "${${variable_to_set}_ORIG}" CACHE INTERNAL "")
    else()
      message("-- Not wrapping ${op_name} because "
        "${variable_to_set}=`${variable_to_set}` is not set."
        " To enable statistics set ${variable_to_set}.")
    endif()
  endif()
endfunction()
#
# NOTE: Above, if this is not the first configure (and
# BUILD_STATS_COMPLETED_FIRST_CONFIG is unset) then we don't want to do
# anything different with the build stats wrappers.  For example, we don't
# want CMAKE_CXX_FLAGS to be empty on the first configure when this function
# is called and have CMake to find the C++ compiler later in the first
# configure and then on the reconfigure have a build stats wrapper generated
# for the C++ compiler.  If this happened, then the C++ code would build with
# the raw C++ compiler after the first configure but after the second and
# subsequent (re)configures would (re)build the code with the build-stats
# wrapped compiler.  It seems like a bad idea to have the code build
# differently on a reconfigure even if the user does not do anything other
# than trigger a reconfigure (e.g. by touching a CMakeLists.txt file or adding
# a new source file).


# Helper macro to shorten above function some
#
# Sets ${variable_to_set} if ${op_lc} is found.
#
macro(generate_build_stats_wrapper_for_op_find_op_lc)
  # there's an issue here - if CMAKE_FOO is unset (whatever `variable_to_set` is)
  # we need a to know the command - but CMake hasn't chosen one yet...
  if( BUILD_STATS_ALLOW_FIND
    AND (
      ("${${variable_to_set}}" STREQUAL "")
      OR
      (NOT ${variable_to_set})
      )
    )
    message("-- " "${variable_to_set} is not set, but a wrapper has been requested. Asking CMake to find ${op_lc}")
    find_program(${variable_to_set} "${op_lc}")
    print_var(${variable_to_set})
  endif()
endmacro()


# Remove the build stats file on each configure if asked to do so.
#
function(remove_build_stats_file_on_configure)

  advanced_set(${PROJECT_NAME}_REMOVE_BUILD_STATS_ON_CONFIGURE OFF
    ${${PROJECT_NAME}_REMOVE_BUILD_STATS_ON_CONFIGURE_DEFAULT} CACHE BOOL
    "If set to 'ON', then the build_stats.csv file will be removed on each configure." )

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


