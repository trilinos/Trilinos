# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

################################################################################
#
# Module containing functions for getting Git repo information
#
################################################################################


# @FUNCTION: tribits_git_repo_sha1()
#
# Get the Git repo version SHA1.
#
# Usage::
#
#   tribits_git_repo_sha1(<gitRepoDir> <gitRepoSha1Out>
#     FAILURE_MESSAGE_OUT <failureMsgOut>)
#
# If ``<failureMsgOut>`` is non-empty on input, then the variable
# ``<failureMsgOut>`` will be set to a non-empty value on output with the
# error message if an error occurs and the SHA1 for ``<gitRepoDir>`` can not
# be computed.  If ``<failureMsgOut>`` is empty in input, and there is an
# error, a ``FATAL_ERROR`` will be raised.
#
# NOTE: To get the SHA1 for the Git repo ``<gitRepoDir>``, it must contain the
# ``.git/`` directory and the ``git`` executable var ``GIT_EXECUTABLE`` must
# be set.  Otherwise, it is an error.
#
function(tribits_git_repo_sha1  gitRepoDir  gitRepoSha1Out)

  cmake_parse_arguments( PARSE_ARGV 2
    PARSE "" "" # prefix, options, one_value_keywords
    "FAILURE_MESSAGE_OUT" # multi_value_keywords
    )
  tribits_check_for_unparsed_arguments()
  tribits_assert_parse_arg_zero_or_one_value(PARSE  FAILURE_MESSAGE_OUT)

  set(failureMsg "")

  if (NOT GIT_EXECUTABLE)
    set(failureMsg "ERROR: The program '${GIT_NAME}' could not be found!")
  elseif (NOT IS_DIRECTORY "${gitRepoDir}/.git")
    set(failureMsg "ERROR: The directory ${gitRepoDir}/.git does not exist!")
  endif()

  set(gitRepoSha1 "")
  if (failureMsg STREQUAL "")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} log -1 "--pretty=format:%H"
      WORKING_DIRECTORY  ${gitRepoDir}
      RESULT_VARIABLE  gitCmndRtn OUTPUT_VARIABLE  gitCmndOutput
      OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_STRIP_TRAILING_WHITESPACE
      )

    if (NOT gitCmndRtn STREQUAL 0)
      set(failureMsg "ERROR: ${GIT_EXECUTABLE} command returned ${gitCmndRtn}!=0 for repo ${gitRepoDir}!")
    else()
      set(gitRepoSha1 "${gitCmndOutput}")
    endif()
  endif()

  if (NOT  failureMsg  STREQUAL "" AND  PARSE_FAILURE_MESSAGE_OUT  STREQUAL "")
    message(FATAL_ERROR "${failureMsg}")
  elseif (PARSE_FAILURE_MESSAGE_OUT)
    set(${PARSE_FAILURE_MESSAGE_OUT} "${failureMsg}" PARENT_SCOPE)
  endif()
  set(${gitRepoSha1Out} "${gitRepoSha1}" PARENT_SCOPE)

endfunction()


# @FUNCTION: tribits_generate_single_repo_version_string()
#
# Get the formatted string containing the current git repo version.
#
# Usage:
#
#   tribits_generate_single_repo_version_string(<gitRepoDir>
#     <repoVersionStringOut> [INCLUDE_PARENT_COMMITS ON|OFF])
#
# If the optional argument ``INCLUDE_PARENT_COMMITS <val>`` is passed,
# then the head commit's parent(s) info will be be included in
# the repo version output string formatted.
#
function(tribits_generate_single_repo_version_string  gitRepoDir
   repoVersionStringOut
  )

  cmake_parse_arguments( PARSE_ARGV 2
    PARSE "" # prefix, optional
    "INCLUDE_COMMIT_PARENTS" "" # one_value_keywords, multi_value_keyword
    )
  tribits_check_for_unparsed_arguments()
  tribits_assert_parse_arg_zero_or_one_value(PARSE INCLUDE_COMMIT_PARENTS)

  tribits_assert_git_executable()

  # A) Get HEAD commit's info
  
  tribits_generate_commit_info_string(${gitRepoDir} HEAD commitInfoString)

  set(outStringBuilder ${commitInfoString})

  # B) Get all of HEAD commit's parents into a list

  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 "--pretty=format:%p" HEAD
    WORKING_DIRECTORY ${gitRepoDir}
    RESULT_VARIABLE gitCmndRtn OUTPUT_VARIABLE gitCmndOutput
    OUTPUT_STRIP_TRAILING_WHITESPACE  ERROR_STRIP_TRAILING_WHITESPACE
  )

  if (NOT gitCmndRtn STREQUAL 0)
    message(FATAL_ERROR "ERROR, ${GIT_EXECUTABLE} command returned ${gitCmndRtn}!=0"
      " with output '${gitCmndOutput}' for sha1 ${gitHeadSha1} of repo ${gitRepoDir}!")
    set(headParentList "Error, could not get commit's parents!")
  else()
    string(REGEX REPLACE " +" ";" headParentList "${gitCmndOutput}")
  endif()

  list(LENGTH headParentList headNumParents)

  # C) Get each parent's commit info and format the output

  if (PARSE_INCLUDE_COMMIT_PARENTS)

    set(parentIdx 1) # Parent commit indexes are 1-based by git

    foreach(parentSha1 IN LISTS headParentList)

      # C.1) Get parent commit info string

      tribits_generate_commit_info_string(
        ${gitRepoDir} ${parentSha1}
        commitInfoString)

      # C.2) Format parent string to be pretty in config output

      string(APPEND outStringBuilder
        "\n    *** Parent ${parentIdx}:")
      string(REPLACE "\n" "\n    "
        commitInfoString "${commitInfoString}")
      string(APPEND outStringBuilder "\n    ${commitInfoString}" )

      math(EXPR parentIdx "${parentIdx}+1")

    endforeach()

  endif()

  set(${repoVersionStringOut} "${outStringBuilder}" PARENT_SCOPE)

endfunction()


# @FUNCTION: tribits_generate_commit_info_string()
#
# Get the formatted commit info containing commit's SHA1,
# author, date, email, and 80 character summary.
#
# Usage:
#   tribits_generate_commit_info_string(<gitRepoDir> <gitCommitSha1>
#     commitInfoStringOut)
#
# NOTE: Below, it is fine if ${maxSummaryLen} > len(${gitCmndOutput}) as
# string(SUBSTRING ...) will just shorten this to the length of the string.
#
function(tribits_generate_commit_info_string gitRepoDir gitCommitSha1
   commitInfoStringOut
  )

  # A) Get commit hash, author, date, and email

  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 "--pretty=format:%h [%ad] <%ae>" ${gitCommitSha1}
    WORKING_DIRECTORY ${gitRepoDir}
    RESULT_VARIABLE gitCmndRtn OUTPUT_VARIABLE gitCmndOutput
    OUTPUT_STRIP_TRAILING_WHITESPACE  ERROR_STRIP_TRAILING_WHITESPACE
    )

  if (NOT gitCmndRtn STREQUAL 0)
    message(FATAL_ERROR "ERROR, ${GIT_EXECUTABLE} command returned ${gitCmndRtn}!=0"
      " with output '${gitCmndOutput}' for sha1 ${gitCommitSha1} of repo ${gitRepoDir}!")
    set(gitVersionLine "Error, could not get version info!")
  else()
    set(gitVersionLine "${gitCmndOutput}")
  endif()

  # B) Get the first 80 chars of the summary message for more info

  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 "--pretty=format:%s" ${gitCommitSha1}
    WORKING_DIRECTORY ${gitRepoDir}
    RESULT_VARIABLE gitCmndRtn OUTPUT_VARIABLE gitCmndOutput
    OUTPUT_STRIP_TRAILING_WHITESPACE  ERROR_STRIP_TRAILING_WHITESPACE
    )

  if (NOT gitCmndRtn STREQUAL 0)
    message(FATAL_ERROR "ERROR, ${GIT_EXECUTABLE} command returned ${gitCmndRtn}!=0"
      " with output '${gitCmndOutput}' for sha1 ${gitCommitSha1} of repo ${gitRepoDir}!")
    set(gitSummaryStr "Error, could not get version summary!")
  else()
    set(maxSummaryLen 80)
    string(SUBSTRING "${gitCmndOutput}" 0 ${maxSummaryLen} gitSummaryStr)
  endif()

  set(${commitInfoStringOut}
    "${gitVersionLine}\n${gitSummaryStr}" PARENT_SCOPE)

endfunction()


function(tribits_assert_git_executable)
  if (NOT GIT_EXECUTABLE)
    message(SEND_ERROR "ERROR, the program '${GIT_NAME}' could not be found!"
      "  We can not generate the repo version file!")
  endif()
endfunction()
