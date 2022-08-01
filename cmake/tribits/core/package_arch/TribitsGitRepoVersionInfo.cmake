# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
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
      RESULT_VARIABLE  gitCmndRtn  OUTPUT_VARIABLE  gitCmndOutput
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


# Run the git log command to get the version info for a git repo
#
function(tribits_generate_single_repo_version_string  gitRepoDir
   repoVersionStringOut
  )

  tribits_assert_git_executable()

  # A) Get the basic version info.

  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 "--pretty=format:%h [%ad] <%ae>"
    WORKING_DIRECTORY ${gitRepoDir}
    RESULT_VARIABLE gitCmndRtn
    OUTPUT_VARIABLE gitCmndOutput
    )

  if (NOT gitCmndRtn STREQUAL 0)
    message(FATAL_ERROR "ERROR, ${GIT_EXECUTABLE} command returned ${gitCmndRtn}!=0"
      " for repo ${gitRepoDir}!")
    set(gitVersionLine "Error, could not get version info!")
  else()
    set(gitVersionLine "${gitCmndOutput}")
  endif()

  # B) Get the first 80 chars of the summary message for more info

  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:%s
    WORKING_DIRECTORY ${gitRepoDir}
    RESULT_VARIABLE gitCmndRtn
    OUTPUT_VARIABLE gitCmndOutput
    )

  if (NOT gitCmndRtn STREQUAL 0)
    message(FATAL_ERROR "ERROR, ${GIT_EXECUTABLE} command returned ${gitCmndRtn}!=0"
      " for extra repo ${gitRepoDir}!")
    set(gitSummaryStr "Error, could not get version summary!")
  else()
    set(maxSummaryLen 80)
    string(SUBSTRING "${gitCmndOutput}" 0 ${maxSummaryLen} gitSummaryStr)
  endif()

  set(${repoVersionStringOut}
    "${gitVersionLine}\n${gitSummaryStr}" PARENT_SCOPE)

endfunction()
# NOTE: Above, it is fine if ${maxSummaryLen} > len(${gitCmndOutput}) as
# string(SUBSTRING ...) will just shorten this to the length of the string.


function(tribits_assert_git_executable)
  if (NOT GIT_EXECUTABLE)
    message(SEND_ERROR "ERROR, the program '${GIT_NAME}' could not be found!"
      "  We can not generate the repo version file!")
  endif()
endfunction()
