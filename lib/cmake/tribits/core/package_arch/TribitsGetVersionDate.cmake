include("${${PROJECT_NAME}_TRIBITS_DIR}/core/utils/MessageWrapper.cmake")
include("${${PROJECT_NAME}_TRIBITS_DIR}/core/utils/TribitsStripQuotesFromStr.cmake")


# @FUNCTION: tribits_get_raw_git_commit_utc_time()
#
# Get the git commit date of a repo at a given commit in UTC in the format
# "2019-03-22 15:34:43 +0000"
#
# Usage::
#
#   tribits_get_raw_git_commit_utc_time(<repo_base_dir> <commit_Ref>
#     <git_commit_utc_time_out> )
#
# This requires find_package(Git) to have been called before calling this
# function and it requires a vesrion of git of 2.10.0 or greater (because
# ``git log`` must support the argument ``--date=iso-local``).
#
function(tribits_get_raw_git_commit_utc_time  repo_base_dir  commit_ref
  git_commit_utc_time_out
  )
  if ("${GIT_EXECUTABLE}" STREQUAL "")
    message_wrapper(FATAL_ERROR "Error, GIT_EXECUTABLE not set!")
  endif()
  if ("${GIT_VERSION_STRING}" STREQUAL "")
    message_wrapper(FATAL_ERROR "Error, GIT_VERSION_STRING not set!")
  endif()
  if (GIT_VERSION_STRING  VERSION_LESS  "2.10.0")
    message_wrapper(FATAL_ERROR
      "Error, GIT_VERSION_STRING=${GIT_VERSION_STRING} < 2.10.0!")
  endif()
  if (NOT TRIBITS_GET_RAW_GIT_COMMIT_UTC_TIME_UNIT_TEST_MODE)
    set(OLD_ENV_TZ "$ENV{TZ}")
    set(ENV{TZ} GMT)
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" log
        "--format=%cd" --date=iso-local -1 ${commit_ref}
      WORKING_DIRECTORY "${repo_base_dir}"
      OUTPUT_VARIABLE  GIT_CMND_OUTPUT
      ERROR_VARIABLE  GIT_CMND_OUTPUT
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE  GIT_CMD_RTN
      TIMEOUT 10 #seconds
    )
    set(ENV{TZ} "${OLD_ENV_TZ}")
    #print_var(ENV{TZ})
    if (NOT GIT_CMD_RTN STREQUAL "0")
      message(FATAL_ERROR
        "ERROR: GIT_CMD_RTN=${GIT_CMD_RTN} != 0!\n"
        "Error Message: ${GIT_CMND_OUTPUT}" )
    endif()
  endif()
  # print_var(GIT_CMND_OUTPUT)
  tribits_strip_quotes_from_str("${GIT_CMND_OUTPUT}" git_commit_no_quotes)
  # ToDo: Assert that the date offset is "+0000" or error out!
  set(${git_commit_utc_time_out} "${git_commit_no_quotes}" PARENT_SCOPE)
endfunction()


# @FUNCTION: tribits_get_version_date_from_raw_git_commit_utc_time()
#
# Takes input of the form "YYYY-MM-DD hh:mm:ss +0000" from the git command::
#
#   git log --format="%cd" --date=iso-local -1 <ref>
# 
# and returns the string integer YYYYMMDDhh.
#
# Usage::
#
#   tribits_get_version_date_from_raw_git_commit_utc_time(
#     ""YYYY-MM-DD hh:mm:ss +0000"  <version_date_var> )
#
# This returns a 10-digit integer ``YYYYMMDDhh`` that should fit in a 32-bit
# integer with a max value of ``2^32 / 2 - 1`` = ``2147483647`` and therefore
# should be good until the last hour of of the last day of the last month of
# the year 2147 (i.e. `2147 12 31 23` = `2147123123`).
#
function(tribits_get_version_date_from_raw_git_commit_utc_time
  git_raw_commit_time_utc  version_date_out
  )
  # Split by spaces first " "
  string(REPLACE " " ";"  git_raw_commit_time_utc_space_array
    "${git_raw_commit_time_utc}")
  #print_var(git_raw_commit_time_utc_space_array)
  list(GET git_raw_commit_time_utc_space_array 0 YYYY_MM_DD) # YYYY-MM-DD
  list(GET git_raw_commit_time_utc_space_array 1 hh_mm_ss)   # hh:mm:ss
  list(GET git_raw_commit_time_utc_space_array 2 utc_offset) # +0000
  #print_var(YYYY_MM_DD)
  #print_var(hh_mm_ss)
  #print_var(utc_offset)
  if (NOT utc_offset STREQUAL "+0000")
    message_wrapper(FATAL_ERROR "ERROR, '${git_raw_commit_time_utc}' is NOT"
      " in UTC which would have offset '+0000'!")
  endif()
  # Split YYYY-MM-DD into its components
  string(REPLACE "-" ";" YYYY_MM_DD_array "${YYYY_MM_DD}")
  list(GET YYYY_MM_DD_array 0 YYYY)
  list(GET YYYY_MM_DD_array 1 MM)
  list(GET YYYY_MM_DD_array 2 DD)
  # Split hh:mm:ss into its components
  string(REPLACE ":" ";" hh_mm_ss_array "${hh_mm_ss}")
  list(GET hh_mm_ss_array 0 hh)
  # Form the full YYYYMMDDhhmm integer and return
  set(${version_date_out} "${YYYY}${MM}${DD}${hh}" PARENT_SCOPE)
endfunction()
