#
# cmake -P script to do and update of the base git repo.
#
# Usage:
#
#  cmake [arguments] -P tribits_ctest_update_commands.cmake
#
# where the -D<var>=<value> arguments are shown below.
#
# The list of commands in this script completely clean out a git repo and
# create a local branch with name ${BRANCH} tracking a remote branch ${BRANCH}
# in the remote repo ${REMOTE_NAME}.  This is robust no matter what the
# current state of the local git repo.  The built-in git commands used by
# ctest_update() are not robust to some use cases where these commands are.
# For example, the commands are robust to the situation where the local repo
# may be tracking a remote branch that may have been deleted in the remote
# repo.  The default commands used in ctest_update() (at least as of CMake
# 3.12) crash in that case.
#

cmake_minimum_required(VERSION 3.23.0 FATAL_ERROR)

message("\ncmake -P tribits_ctest_update_commands.cmake:")
message("-- GIT_EXE=${GIT_EXE}")
message("-- REMOTE_NAME=${REMOTE_NAME}")
message("-- BRANCH=${BRANCH}")
message("-- UNIT_TEST_MODE=${UNIT_TEST_MODE}")

set(OVERALL_SUCCESS TRUE)
set(ERROR_CODE 0)

macro(execute_process_wrapper)
  if (UNIT_TEST_MODE)
    message("execute_process(${ARGN})")
  else()
    execute_process(${ARGN} RESULT_VARIABLE RTN_CODE)
    message("RTN_CODE: ${RTN_CODE}")
    if (NOT "${RTN_CODE}" STREQUAL "0")
      set(OVERALL_SUCCESS FALSE)
      set(ERROR_CODE ${RTN_CODE})
    endif()
  endif()
endmacro()

macro(run_command)
  string(REPLACE ";" " " CMND_STR "${ARGN}")
  message("\nRunning: ${CMND_STR}")
  execute_process_wrapper(COMMAND ${ARGN})
endmacro()

run_command(
  "${GIT_EXE}" fetch ${REMOTE_NAME} )

run_command(
  "${GIT_EXE}" clean -fdx )

run_command(
  "${GIT_EXE}" reset --hard HEAD )

if (BRANCH)
  run_command(
    "${GIT_EXE}" checkout -B ${BRANCH} --track ${REMOTE_NAME}/${BRANCH} )
else()
  run_command(
    "${GIT_EXE}" reset --hard @{u} )
endif()

if (OVERALL_SUCCESS)
  message("\nGit Update PASSED!")
else()
  message(FATAL_ERROR "Git Update FAILED!")
endif()

# NOTE: Above, you have to use separate execute_process() commands for each
# git command or you get git errors complaining about git commands running on
# top of each other.  The execute_process() implementation must not ensure
# that one command is completely finished before the next one starts.
