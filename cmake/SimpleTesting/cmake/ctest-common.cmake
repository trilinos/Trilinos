include_guard()
message("+--------------------------------------+")
message("| ctest-common.cmake START             |")
message("+--------------------------------------+")

cmake_minimum_required(VERSION 3.17.0 FATAL_ERROR)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
include(${CMAKE_CURRENT_LIST_DIR}/ctest-functions.cmake)


# -----------------------------------------------------------
# -- Required Parameters
# -----------------------------------------------------------
if(DEFINED build_name)
    set(CTEST_BUILD_NAME ${build_name})
else()
    message(FATAL_ERROR "No build_name supplied. This is necessary to run and submit results.")
endif()

# `subprojects_file` is a -D parameter that is passed in from Trilinos
# See: https://github.com/trilinos/Trilinos/blob/master/cmake/std/trilinosprhelpers/TrilinosPRConfigurationStandard.py
include( ${subprojects_file} OPTIONAL RESULT_VARIABLE subproject_include_success )
if( NOT subproject_include_success )
    message(FATAL_ERROR "Subproject label file '${subprojects_file}' not found.")
else()
    message("Subproject label file '${subprojects_file}' loaded.")
endif()
list(LENGTH CTEST_LABELS_FOR_SUBPROJECTS subproject_count)

# configure_script is passed into CTest via -Dconfigure_script, points to
# Trilinos/cmake/std/<config_script>
get_filename_component(configure_file ${configure_script} ABSOLUTE)

# package_enables is passed into CTest via -Dpackage_enables, points to the
# generated "PackageEnables.cmake file"
get_filename_component(package_enables_file ${package_enables} ABSOLUTE)


# -----------------------------------------------------------
# -- Optional Parameters
# -----------------------------------------------------------
if( NOT DEFINED ctest_submit_retry_count )
    set( ctest_submit_retry_count 5 )
endif()

if( NOT DEFINED ctest_submit_retry_delay )
    set( ctest_submit_retry_delay 3 )
endif()

if( NOT DEFINED dashboard_model )
    set( dashboard_model Experimental )
endif()

if( NOT DEFINED dashboard_track )
    set( dashboard_track Experimental )
endif()

if( NOT DEFINED skip_by_parts_submit )
    set( skip_by_parts_submit ON )
endif()

if( NOT DEFINED skip_clean_build_dir )
    set( skip_clean_build_dir ON )
endif()

if( NOT DEFINED skip_single_submit )
    set( skip_single_submit ON )
endif()

if( NOT DEFINED skip_update_step )
    set( skip_update_step OFF )
endif()

if( NOT DEFINED skip_upload_config_files )
    set( skip_upload_config_files OFF )
endif()


# -----------------------------------------------------------
# -- Miscellaneous Settings
# -----------------------------------------------------------

# Set the output to English
set($ENV{LC_MESSAGES} "en_EN")

# Set the machine_name
cmake_host_system_information(RESULT machine_name QUERY HOSTNAME)
message(">>> machine_name : ${machine_name}")


# -----------------------------------------------------------
# -- Set Git command
# -----------------------------------------------------------
find_program(CTEST_GIT_COMMAND NAMES git)
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

# -----------------------------------------------------------
# -- build specific
# -----------------------------------------------------------

# -- Set `CTEST_SOURCE_DIRECTORY` which would be the directory that contains the Trilinos repository.
#    Note-1: CTEST_SOURCE_DIRECTORY is the path to the **this** cmake file.
#    Note-2: use of `get_filename_component()` is just getting the path up 1
#            level.
if(DEFINED source_dir)
    set(CTEST_SOURCE_DIRECTORY ${source_dir})
else()
    message(FATAL_ERROR "Missing required parameter: `source_dir`")
endif()
message(">>> CTEST_SOURCE_DIRECTORY = ${CTEST_SOURCE_DIRECTORY}")


## -- build_root
if( NOT DEFINED build_root )
    get_filename_component(bin_dir ${CTEST_SOURCE_DIRECTORY} DIRECTORY)
    set(build_root "${bin_dir}/nightly_testing")
endif()


# -- build_dir
if(DEFINED build_dir)
    message(">>> set CTEST_BINARY_DIRECTORY = ${build_dir} (A)")
    set(CTEST_BINARY_DIRECTORY "${build_dir}")
else()
    message(">>> set CTEST_BINARY_DIRECTORY = ${build_root}/${CTEST_BUILD_NAME} (B)")
    set(CTEST_BINARY_DIRECTORY "${build_root}/${CTEST_BUILD_NAME}")
endif()
message(">>> CTEST_BINARY_DIRECTORY = ${CTEST_BINARY_DIRECTORY}")


# -----------------------------------------------------------
# -- Configure parallelism
# -----------------------------------------------------------
if( NOT DEFINED PARALLEL_LEVEL )
    cmake_host_system_information(RESULT PARALLEL_LEVEL QUERY NUMBER_OF_PHYSICAL_CORES)
endif()

if( NOT DEFINED TEST_PARALLEL_LEVEL )
    set(TEST_PARALLEL_LEVEL "${PARALLEL_LEVEL}")
endif()

# These are command line options to Ninja
# See: https://manpages.debian.org/testing/ninja-build/ninja.1.en.html
# -j${PARALLEL_LEVEL} - now many concurrent build jobs to allow
# -k NNN - keep going until NNN jobs fail. 0 == infinity
# -l N   - might be an interesting option to reduce # of jobs based on LOAD AVERAGE
#          but load average is dependent on the # of cores on a system.
set(CTEST_BUILD_FLAGS "-j${PARALLEL_LEVEL} -k 0")




# -----------------------------------------------------------
# -- Dependent Variables
# -----------------------------------------------------------

# * REQUIRES `CTEST_BINARY_DIRECTORY` to be set.
if(NOT skip_upload_config_files)
    set(configure_command_file ${CTEST_BINARY_DIRECTORY}/configure_command.txt)
    set(genconfig_build_name_file ${CTEST_BINARY_DIRECTORY}/genconfig_build_name.txt)
endif()


# -----------------------------------------------------------
# -- CTest Settings
# -----------------------------------------------------------

# Set CTEST_SITE to the name of the system.
cmake_host_system_information(RESULT HOSTNAME QUERY HOSTNAME)

set(CTEST_SITE "${HOSTNAME}")

# See: https://cmake.org/cmake/help/latest/command/site_name.html#command:site_name
site_name(${CTEST_SITE})

# This is required to make build information detailed enough to
# enable assignment to a specific subproject.
# This must be set for both cmake and ctest to work.
set(CTEST_USE_LAUNCHERS ON)
set(ENV{CTEST_USE_LAUNCHERS_DEFAULT} 1)


# -----------------------------------------------------------
# -- Print out settings
# -----------------------------------------------------------
print_options_list()


# -----------------------------------------------------------
# -- Set up CDash
# -----------------------------------------------------------
# include(${CMAKE_CURRENT_LIST_DIR}/ctest-cdash-setup.cmake)

message("+--------------------------------------+")
message("| ctest-common.cmake FINISH            |")
message("+--------------------------------------+")
