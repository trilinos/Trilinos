set(CTEST_SITE "$ENV{CTEST_SITE}")
if("${CTEST_SITE}" STREQUAL "")
  site_name(CTEST_SITE)
endif()
if("${CTEST_SITE}" STREQUAL "")
  if(WIN32)
    string(TOLOWER "$ENV{COMPUTERNAME}" CTEST_SITE)
  else()
    execute_process(COMMAND uname -n
      OUTPUT_VARIABLE CTEST_SITE
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
  endif()
endif()

set(CTEST_BUILD_NAME "$ENV{CTEST_BUILD_NAME}")
if("${CTEST_BUILD_NAME}" STREQUAL "")
  if(WIN32)
    set(HOST_TYPE $ENV{OS})
  else()
    execute_process(COMMAND uname
      OUTPUT_VARIABLE HOST_TYPE
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
  endif()
  set(CTEST_BUILD_NAME "${HOST_TYPE}-TDD-${CTEST_SITE}")
endif()

set(CTEST_CMAKE_GENERATOR "$ENV{CTEST_CMAKE_GENERATOR}")
if("${CTEST_CMAKE_GENERATOR}" STREQUAL "")
  if(WIN32)
    set(CTEST_CMAKE_GENERATOR "NMake Makefiles")
  else()
    set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
  endif()
endif()

set(CTEST_TEST_TIMEOUT "$ENV{CTEST_TEST_TIMEOUT}")
if("${CTEST_TEST_TIMEOUT}" STREQUAL "")
  set(CTEST_TEST_TIMEOUT 7200)
endif()

get_filename_component(CTEST_SOURCE_DIRECTORY
  "${CTEST_SCRIPT_DIRECTORY}" ABSOLUTE)

get_filename_component(CTEST_BINARY_DIRECTORY
  "${CTEST_SCRIPT_DIRECTORY}/../../../../TDD_BUILD" ABSOLUTE)

get_filename_component(CTEST_NOTES_FILES
  "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" ABSOLUTE)
if(NOT "$ENV{TDD_CRON_DRIVER_LOGFILE}" STREQUAL "")
  set(CTEST_NOTES_FILES ${CTEST_NOTES_FILES} "$ENV{TDD_CRON_DRIVER_LOGFILE}")
endif()
if(NOT "$ENV{TDD_CRON_DRIVER_SCRIPT}" STREQUAL "")
  set(CTEST_NOTES_FILES ${CTEST_NOTES_FILES} "$ENV{TDD_CRON_DRIVER_SCRIPT}")
endif()

set(parallel_level "$ENV{TDD_PARALLEL_LEVEL}")
if("${parallel_level}" STREQUAL "")
  set(parallel_level 1)
endif()

set(git_exe "$ENV{TDD_GIT_EXE}")
if("${git_exe}" STREQUAL "")
  find_program(git_exe NAMES git.cmd eg git)
endif()
if(git_exe)
  set(CTEST_UPDATE_COMMAND "${git_exe}")
endif()

ctest_empty_binary_directory("${CTEST_BINARY_DIRECTORY}")
ctest_start("Experimental")
ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}")
ctest_configure(BUILD "${CTEST_BINARY_DIRECTORY}")
ctest_read_custom_files(BUILD "${CTEST_BINARY_DIRECTORY}")
ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}")
ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" PARALLEL_LEVEL ${parallel_level})
ctest_submit()
