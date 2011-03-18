
MESSAGE(
 "\n***"
 "\n*** TrilinosDriverDashboard.cmake"
 "\n***\n"
 )

# Get the base diretory for the Trilinos source.  We only assume that the
# CTest script that is being called is under Trilinos/cmake.
STRING(REGEX MATCH "(.+/Trilinos)/cmake" TRILINOS_CMAKE_DIR
  "${CTEST_SCRIPT_DIRECTORY}" )
IF("${TRILINOS_CMAKE_DIR}" STREQUAL "")
  STRING(REGEX MATCH "(.+)/cmake" TRILINOS_CMAKE_DIR
    "${CTEST_SCRIPT_DIRECTORY}" )
ENDIF()
  
MESSAGE("TRILINOS_CMAKE_DIR = ${TRILINOS_CMAKE_DIR}")

SET( CMAKE_MODULE_PATH
   "${TRILINOS_CMAKE_DIR}/utils"
   )

#MESSAGE("CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")

INCLUDE(PrintVar)
INCLUDE(SetDefaultAndFromEnv)

#
# A) Set up the environment get options
#

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
  
# Submit the results to the dashboard or not
SET_DEFAULT_AND_FROM_ENV( TDD_DO_SUBMIT TRUE )

# Dashboard model : Nightly, Experimental, Continuous
SET_DEFAULT_AND_FROM_ENV( TDD_CTEST_TEST_TYPE Experimental )

# set this to ON if you need to test something before committing.
SET_DEFAULT_AND_FROM_ENV( TDD_IN_TESTING_MODE OFF )

get_filename_component(CTEST_SOURCE_DIRECTORY
  "${CTEST_SCRIPT_DIRECTORY}" ABSOLUTE)

get_filename_component(CTEST_UPDATE_DIRECTORY
  "${CTEST_SCRIPT_DIRECTORY}/../../.." ABSOLUTE)

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
  set(git_exe "git_exe-NOTFOUND")
  find_program(git_exe NAMES git.cmd eg git)
endif()
if(git_exe)
  set(CTEST_UPDATE_TYPE "git")
  set(CTEST_UPDATE_COMMAND "${git_exe}")
endif()

#
# Run the outer dashboard
#

message("\nA) Empty out ${CTEST_BINARY_DIRECTORY} ...")
ctest_empty_binary_directory("${CTEST_BINARY_DIRECTORY}")

ctest_start("${TDD_CTEST_TEST_TYPE}")

message("\nB) Update ${CTEST_UPDATE_DIRECTORY} ...")
message("      CTEST_UPDATE_COMMAND='${CTEST_UPDATE_COMMAND}'")
message("      CTEST_UPDATE_TYPE='${CTEST_UPDATE_TYPE}'")
if(NOT TDD_IN_TESTING_MODE)
  ctest_update(SOURCE "${CTEST_UPDATE_DIRECTORY}")
else()
  message("In testing mode no update will be performed.")
endif()

message("\nC) Configure ${CTEST_BINARY_DIRECTORY} ...")
ctest_configure(BUILD "${CTEST_BINARY_DIRECTORY}")

ctest_read_custom_files(BUILD "${CTEST_BINARY_DIRECTORY}")

message("\nD) Build ${CTEST_BINARY_DIRECTORY} ...")
ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND)

message("\nE) Submitting update configure notes build ...")
if (TDD_DO_SUBMIT)
  if(NOT "$ENV{TDD_CTEST_DROP_SITE}" STREQUAL "")
    set(CTEST_DROP_SITE "$ENV{TDD_CTEST_DROP_SITE}")
  endif()
  if(NOT "$ENV{TDD_CTEST_DROP_LOCATION}" STREQUAL "")
    set(CTEST_DROP_LOCATION "$ENV{TDD_CTEST_DROP_LOCATION}")
  endif()
  ctest_submit(PARTS update configure notes build)
else()
  message("\nSkipping submit!")
endif()

message("\nF) Run tests (which run all everything really): PARALLEL_LEVEL ${parallel_level} from ${CTEST_BINARY_DIRECTORY} ...")
ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" PARALLEL_LEVEL ${parallel_level})

message("\nG) Submitting Test ...")
if (TDD_DO_SUBMIT)
  ctest_submit(PARTS Test)
else()
  message("\nSkipping submit!")
endif()

MESSAGE(
 "\n*** Finished TrilinosDriverDashboard.cmake ***\n"
 )
