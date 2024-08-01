# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


include(PrintVar)
include(TribitsParseArgumentsHelpers)
include(TribitsReadTagFile)


################################################################################
#
# This module contains functions for constructing CDash URLs to build and test
# results from inside of a CTest -S script.
#
################################################################################


# @FUNCTION: tribits_get_cdash_results_string_and_write_to_file()
#
# Calls `tribits_get_cdash_results_urls_string()`_ and then writes the CDash
# URLs to a file.
#
# Usage::
#
#   tribits_get_cdash_results_string_and_write_to_file(
#     [CDASH_RESULTS_STRING_OUT <cdashResultsStringOut>]
#     [CDASH_RESULTS_FILE_OUT <cdashResultsFileOut>]
#     )
#
function(tribits_get_cdash_results_string_and_write_to_file)
  # Parse args
  cmake_parse_arguments( PARSE_ARGV 0
    PARSE "" "" # prefix, options, one_value_keywords
    "CDASH_RESULTS_STRING_OUT;CDASH_RESULTS_FILE_OUT" # multi_value_keywords
    )
  tribits_check_for_unparsed_arguments(PARSE)
  tribits_assert_parse_arg_zero_or_one_value(PARSE  CDASH_RESULTS_STRING_OUT
    CDASH_RESULTS_FILE_OUT )
  # Get and set CDash results URL
  tribits_get_cdash_results_urls_string(cdashResultsString)
  if (PARSE_CDASH_RESULTS_STRING_OUT)
    set(${PARSE_CDASH_RESULTS_STRING_OUT} "${cdashResultsString}" PARENT_SCOPE)
  endif()
  if (PARSE_CDASH_RESULTS_FILE_OUT)
    file(WRITE "${PARSE_CDASH_RESULTS_FILE_OUT}" "${cdashResultsString}")
  endif()
endfunction()


# @FUNCTION: tribits_get_cdash_results_urls_string()
#
# Call `tribits_get_cdash_results_urls()`_ and then construct a CDash URLs
# string fit for printing.
#
# Usage::
#
#   tribits_get_cdash_results_urls_string(
#     <cdashResultsUrlsStringOut>
#     )
#
# Construct the build and test URLs on CDash given the site name, buildname,
# and buildstamp (taken from the TAG file) from inside of a running ctest -S
# program and optionally write it to a file as well.
#
function(tribits_get_cdash_results_urls_string  cdashResultsUrlsStringOut)
  tribits_get_cdash_results_urls(
    CDASH_BUILD_URL_OUT  cdashBuildUrl
    CDASH_REVISION_BUILDS_URL_OUT  cdashRevisionBuildsUrl
    CDASH_REVISION_NONPASSING_TESTS_URL_OUT  cdashRevisionNonpassingTestsUrl
    )
  tribits_generate_cdash_results_string_from_urls(
    CDASH_BUILD_URL  "${cdashBuildUrl}"
    CDASH_REVISION_BUILDS_URL  "${cdashRevisionBuildsUrl}"
    CDASH_REVISION_NONPASSING_TESTS_URL  "${cdashRevisionNonpassingTestsUrl}"
    CDASH_RESULTS_STRING_OUT  cdashResultsUrlsString
    )
  set(${cdashResultsUrlsStringOut} "${cdashResultsUrlsString}" PARENT_SCOPE)
endfunction()


# @FUNCTION: tribits_generate_cdash_results_string_from_urls()
#
# Generate the CDash results URL string given the individual URLs.
#
# Usage::
#
#   tribits_generate_cdash_results_string_from_urls(
#     CDASH_BUILD_URL  "<cdashBuildUrl>"
#     [CDASH_REVISION_BUILDS_URL  "<cdashRevisionBuildsUrl>"]
#     [CDASH_REVISION_NONPASSING_TESTS_URL  "<cdashRevisionNonpassingTestsUrl>"]
#     CDASH_RESULTS_STRING_OUT  <cdashResultsUrlsStringOut>
#     )
#
# Takes the URLs returned from `tribits_get_cdash_results_urls()`_ and
# generates a string out of them which is set in the return var
# ``<cdashResultsUrlsStringOut>``.
#
function(tribits_generate_cdash_results_string_from_urls)
  # Parse args
  cmake_parse_arguments(PARSE_ARGV 0
    PARSE "" "" # prefix, options, one_value_keywords
    # multi_value_keywords
    "CDASH_BUILD_URL;CDASH_REVISION_BUILDS_URL;CDASH_REVISION_NONPASSING_TESTS_URL;CDASH_RESULTS_STRING_OUT"
    )
  tribits_check_for_unparsed_arguments()
  tribits_assert_parse_arg_one_value(PARSE  CDASH_BUILD_URL
    CDASH_RESULTS_STRING_OUT)
  tribits_assert_parse_arg_zero_or_one_value(PARSE  CDASH_REVISION_BUILDS_URL
    CDASH_REVISION_NONPASSING_TESTS_URL)
  # Construct CDash results URLs string
  set(cdashResultsString "")
  string(APPEND  cdashResultsString
    "Link to this build's results on CDash:\n"
    "\n"
    "    ${PARSE_CDASH_BUILD_URL}\n")
  if (PARSE_CDASH_REVISION_BUILDS_URL)
    string(APPEND  cdashResultsString
      "\nLink to all builds for this repo version on CDash:\n"
      "\n"
      "    ${PARSE_CDASH_REVISION_BUILDS_URL}\n")
  endif()
  if (PARSE_CDASH_REVISION_NONPASSING_TESTS_URL)
    string(APPEND  cdashResultsString
      "\nLink to all nonpassing tests for all builds for this repo version on CDash:\n"
      "\n"
      "    ${PARSE_CDASH_REVISION_NONPASSING_TESTS_URL}\n")
  endif()
  # Set output
  set(${PARSE_CDASH_RESULTS_STRING_OUT} ${cdashResultsString} PARENT_SCOPE)
endfunction()


# @FUNCTION: tribits_get_cdash_results_urls()
#
# Construct the build and test URLs on CDash given the site name, buildname,
# and buildstamp (taken from the TAG file) from inside of a running ctest -S
# program and optionally write it to a file as well.
#
# Usage::
#
#   tribits_get_cdash_results_urls(
#     CDASH_BUILD_URL_OUT <cdashBuildUrlOut>
#     [CDASH_REVISION_BUILDS_URL_OUT <cdashRevisionBuildsUrlOut>]
#     [CDASH_REVISION_NONPASSING_TESTS_URL_OUT <cdashRevisionNonpassingTestsUrlOut>]
#     )
#
# Here, the CDash URLs are constructed the following CMake variables already
# set in a ``ctest -S`` process:
#
#   * ``CTEST_DROP_SITE``
#   * ``CTEST_DROP_LOCATION`` (``submit.php`` is replaced with ``index.php``)
#   * ``CTEST_PROJECT_NAME``
#   * ``CTEST_SITE``
#   * ``CTEST_BUILD_NAME``
#   * ``CTEST_BINARY_DIRECTORY``
#   * ``CTEST_SOURCE_DIRECTORY``
#
# and other information derived from that.
#
# The buildstamp is read in from the file
# ``${CTEST_BINARY_DIRECTORY}/Testing/TAG``.
#
# If available, the revision SHA1 is obtained from the git repo at
# ``CTEST_SOURCE_DIRECTORY`` if the directory
# ``${CTEST_SOURCE_DIRECTORY}/.git`` exists.  If the base project source
# directory is not a git reposistory, then ``<cdashRevisionBuildsUrlOut>`` and
# ``<cdashRevisionNonpassingTestsUrlOut>``, if requested, will be set to
# empty.
#
# Note that the CDash URLs will have ``https://`` added to the beginning so
# that GitHub Actions and other systems will put in a hyperlink to them.
#
function(tribits_get_cdash_results_urls)
  # Parse args
  cmake_parse_arguments(PARSE_ARGV 0
    PARSE "" "" # prefix, options, one_value_keywords
    # multi_value_keywords
    "CDASH_BUILD_URL_OUT;CDASH_REVISION_BUILDS_URL_OUT;CDASH_REVISION_NONPASSING_TESTS_URL_OUT"
    )
  tribits_check_for_unparsed_arguments(PARSE)
  tribits_assert_parse_arg_one_value(PARSE  CDASH_BUILD_URL_OUT)
  tribits_assert_parse_arg_zero_or_one_value(PARSE  CDASH_REVISION_BUILDS_URL_OUT
    CDASH_REVISION_NONPASSING_TESTS_URL_OUT)
  # Get the info
  tribits_get_cdash_build_url(cdashBuildUrl)
  tribits_git_repo_sha1("${CTEST_SOURCE_DIRECTORY}"  gitRepoSha1
    FAILURE_MESSAGE_OUT  gitRepoSha1FailureMsg)
  if (gitRepoSha1)
    tribits_get_cdash_site_from_drop_site_and_location(
      CTEST_DROP_SITE "${CTEST_DROP_SITE}"
      CTEST_DROP_LOCATION "${CTEST_DROP_LOCATION}"
      CDASH_SITE_URL_OUT  cdashSiteUrl
      )
    tribits_get_cdash_revision_builds_url(
      CDASH_SITE_URL "${cdashSiteUrl}"
      PROJECT_NAME "${CTEST_PROJECT_NAME}"
      GIT_REPO_SHA1  "${gitRepoSha1}"
      CDASH_REVISION_BUILDS_URL_OUT  cdashRevisionBuildsUrl
     )
   tribits_get_cdash_revision_nonpassing_tests_url(
      CDASH_SITE_URL "${cdashSiteUrl}"
      PROJECT_NAME "${CTEST_PROJECT_NAME}"
      GIT_REPO_SHA1  "${gitRepoSha1}"
      CDASH_REVISION_NONPASSING_TESTS_URL_OUT  cdashRevisionNonpassingTestsUrl
      )
  else()
    set(cdashRevisionBuildsUrl "")
    set(cdashRevisionNonpassingTestsUrl "")
  endif()
  # Set the outputs
  set(${PARSE_CDASH_BUILD_URL_OUT} "${cdashBuildUrl}" PARENT_SCOPE)
  set(${PARSE_CDASH_REVISION_BUILDS_URL_OUT} "${cdashRevisionBuildsUrl}"
    PARENT_SCOPE)
  set(${PARSE_CDASH_REVISION_NONPASSING_TESTS_URL_OUT} "${cdashRevisionNonpassingTestsUrl}"
    PARENT_SCOPE)
endfunction()


# @FUNCTION: tribits_get_cdash_revision_builds_url()
#
# Get the CDash URL for all builds with the same repo version SHA1
#
# Usage::
#
#   tribits_get_cdash_revision_builds_url(
#     CDASH_SITE_URL <cdashSiteUrl>
#     PROJECT_NAME <projectName>
#     GIT_REPO_SHA1 <gitRepoSha1>
#     CDASH_REVISION_BUILDS_URL_OUT <cdashRevisionBuildsUrlOut>
#     )
#
function(tribits_get_cdash_revision_builds_url)
  cmake_parse_arguments(PARSE_ARGV 0
    PARSE "" "" # prefix, options, one_value_keywords
    # multi_value_keywords
    "CDASH_SITE_URL;PROJECT_NAME;GIT_REPO_SHA1;CDASH_REVISION_BUILDS_URL_OUT"
    )
  tribits_check_for_unparsed_arguments()
  tribits_assert_parse_arg_one_value(PARSE  CDASH_SITE_URL  PROJECT_NAME
    GIT_REPO_SHA1  CDASH_REVISION_BUILDS_URL_OUT)
  set(${PARSE_CDASH_REVISION_BUILDS_URL_OUT}
    "${PARSE_CDASH_SITE_URL}/index.php?project=${PARSE_PROJECT_NAME}&filtercount=1&showfilters=1&field1=revision&compare1=61&value1=${PARSE_GIT_REPO_SHA1}"
    PARENT_SCOPE)
endfunction()


# @FUNCTION: tribits_get_cdash_revision_nonpassing_tests_url()
#
# Get the CDash URL for all non-passing tests with the same repo version SHA1
#
# Usage::
#
#   tribits_get_cdash_revision_nonpassing_tests_url(
#     CDASH_SITE_URL <cdashSiteUrl>
#     PROJECT_NAME <projectName>
#     GIT_REPO_SHA1 <gitRepoSha1>
#     CDASH_REVISION_NONPASSING_TESTS_URL_OUT <cdashRevisionNonpassingTestsUrlOut>
#     )
#
function(tribits_get_cdash_revision_nonpassing_tests_url)
  cmake_parse_arguments( PARSE_ARGV 0
    PARSE "" "" # prefix, options, one_value_keywords
    # multi_value_keywords
    "CDASH_SITE_URL;PROJECT_NAME;GIT_REPO_SHA1;CDASH_REVISION_NONPASSING_TESTS_URL_OUT"
    )
  tribits_check_for_unparsed_arguments()
  tribits_assert_parse_arg_one_value(PARSE  CDASH_SITE_URL  PROJECT_NAME
    GIT_REPO_SHA1  CDASH_REVISION_NONPASSING_TESTS_URL_OUT)
  set(${PARSE_CDASH_REVISION_NONPASSING_TESTS_URL_OUT}
    "${PARSE_CDASH_SITE_URL}/queryTests.php?project=${PARSE_PROJECT_NAME}&filtercount=2&showfilters=1&filtercombine=and&field1=revision&compare1=61&value1=${PARSE_GIT_REPO_SHA1}&field2=status&compare2=62&value2=passed"
    PARENT_SCOPE)
endfunction()


# @FUNCTION: tribits_get_cdash_build_url()
#
# Construct the build URL on CDash given the site name, buildname, and
# buildstamp (taken from the TAG file) from inside of a running ctest -S
# program.
#
# Usage::
#
#   tribits_get_cdash_build_url(<cdashBuildUrlOut>)
#
# Here, ``<cdashBuildUrlOut>`` returns the CDash Build URL constructed from
# the following CMake variables already set in a ``ctest -S`` process:
#
#   * ``CTEST_DROP_SITE``
#   * ``CTEST_DROP_LOCATION`` (``submit.php`` is replaced with ``index.php``)
#   * ``CTEST_PROJECT_NAME``
#   * ``CTEST_SITE``
#   * ``CTEST_BUILD_NAME``
#
# and the buildstamp read in from the file
# ``${CTEST_BINARY_DIRECTORY}/Testing/TAG``.
#
# Note that ``<cdashBuildUrlOut>`` will have ``https://`` added to the
# beginning of it so that GitHub Actions and other systems will put in a link
# to them.
#
function(tribits_get_cdash_build_url  cdashBuildUrlOut)
  tribits_get_cdash_index_php_from_drop_site_and_location(
    CTEST_DROP_SITE "${CTEST_DROP_SITE}"
    CTEST_DROP_LOCATION "${CTEST_DROP_LOCATION}"
    INDEX_PHP_URL_OUT indexPhpUrl
    )
  tribits_get_cdash_build_url_from_tag_file(
    INDEX_PHP_URL "${indexPhpUrl}"
    PROJECT_NAME "${CTEST_PROJECT_NAME}"
    SITE_NAME "${CTEST_SITE}"
    BUILD_NAME "${CTEST_BUILD_NAME}"
    TAG_FILE "${CTEST_BINARY_DIRECTORY}/Testing/TAG"
    CDASH_BUILD_URL_OUT cdashBuildUrl
    )
  set(${cdashBuildUrlOut} "${cdashBuildUrl}" PARENT_SCOPE)
endfunction()


# @FUNCTION: tribits_print_cdash_url()
#
# Print the URL on CDash where build results can be found.
#
# Usage::
#
#   tribits_print_cdash_url( <msg> <cdashUrl> )
#
function(tribits_print_cdash_url  msg  cdashUrl)
  message("\n${msg}\n")
  message("    ${cdashUrl}\n")
endfunction()


# @FUNCTION: tribits_get_cdash_build_url_from_tag_file()
#
# Create CDash index.php URL from the build parts.
#
# Usage::
#
#   tribits_get_cdash_build_url_from_tag_file(
#     INDEX_PHP_URL <indexPhpUrl>
#     PROJECT_NAME <projectName>
#     SITE_NAME <siteName>
#     BUILD_NAME <buildName>
#     TAG_FILE <tagFile>
#     CDASH_BUILD_URL_OUT <cdashBuildUrlOut>
#     )
#
# Note that spaces are allowed ``<siteName>`` or ``<buildName>`` and those
# will be handled correctly to produce a valid URL.
#
function(tribits_get_cdash_build_url_from_tag_file)
  # Get arguments
  cmake_parse_arguments(
    PREFIX #prefix
    "" #options
    "INDEX_PHP_URL;PROJECT_NAME;SITE_NAME;BUILD_NAME;TAG_FILE;CDASH_BUILD_URL_OUT" #one_value_keywords
    "" #multi_value_keywords
    ${ARGN}
    )
  # Read in the tag file and get the build stamp from that
  tribits_read_ctest_tag_file(${PREFIX_TAG_FILE} buildStartTime cdashGroup
    cdashModel # The model is not used here but we still need to include this arg
    )
  set(buildstamp "${buildStartTime}-${cdashGroup}")
  # Build the URL and return it
  tribits_get_cdash_build_url_from_parts(
    INDEX_PHP_URL "${PREFIX_INDEX_PHP_URL}"
    PROJECT_NAME "${PREFIX_PROJECT_NAME}"
    SITE_NAME "${PREFIX_SITE_NAME}"
    BUILD_NAME "${PREFIX_BUILD_NAME}"
    BUILD_STAMP "${buildstamp}"
    CDASH_BUILD_URL_OUT cdashBuildUrl
    )
  set(${PREFIX_CDASH_BUILD_URL_OUT} "${cdashBuildUrl}" PARENT_SCOPE)
endfunction()


# @FUNCTION: tribits_get_cdash_index_php_from_drop_site_and_location()
#
# Get the CDash index.php URL from the input CTEST_DROP_SITE and
# CTEST_DROP_LOCATION vars used in a ctest -S script.
#
function(tribits_get_cdash_index_php_from_drop_site_and_location)
  cmake_parse_arguments(
    PREFIX #prefix
    "" #options
    "CTEST_DROP_SITE;CTEST_DROP_LOCATION;INDEX_PHP_URL_OUT" #one_value_keywords
    "" #multi_value_keywords
    ${ARGN}
    )
  tribits_get_cdash_site_from_drop_site_and_location(
    CTEST_DROP_SITE  ${PREFIX_CTEST_DROP_SITE}
    CTEST_DROP_LOCATION  ${PREFIX_CTEST_DROP_LOCATION}
    CDASH_SITE_URL_OUT  cdashSiteUrl )
  set(${PREFIX_INDEX_PHP_URL_OUT} "${cdashSiteUrl}/index.php" PARENT_SCOPE)
endfunction()


# @FUNCTION: tribits_get_cdash_site_from_drop_site_and_location()
#
# Get the full CDash site base URL from the input CTEST_DROP_SITE and
# CTEST_DROP_LOCATION vars used in a ctest -S script.
#
# Usage::
#
#   tribits_get_cdash_site_from_drop_site_and_location(
#     CTEST_DROP_SITE <ctestDropSite>
#     CTEST_DROP_LOCATION  <ctestDropLocation>
#     CDASH_SITE_URL_OUT  <cdashSiteUrlOut>
#     )
#
function(tribits_get_cdash_site_from_drop_site_and_location)
  # Parse args
  cmake_parse_arguments(PARSE_ARGV 0
    PREFIX #prefix
    "" #options
    "CTEST_DROP_SITE;CTEST_DROP_LOCATION;CDASH_SITE_URL_OUT" #one_value_keywords
    "" #multi_value_keywords
    )
  tribits_check_for_unparsed_arguments(PREFIX)
  tribits_assert_parse_arg_one_value(PREFIX  CTEST_DROP_SITE)
  tribits_assert_parse_arg_one_value(PREFIX  CTEST_DROP_LOCATION)
  tribits_assert_parse_arg_one_value(PREFIX  CDASH_SITE_URL_OUT)
  # Get the full CDash site from parts
  string(FIND "${PREFIX_CTEST_DROP_LOCATION}" "?" beginningOfQueryString)
  string(SUBSTRING "${PREFIX_CTEST_DROP_LOCATION}" 0 ${beginningOfQueryString} submitPhpPart)
  string(REPLACE "/submit.php" "" endCDashUrl "${submitPhpPart}")
  set(cdashSiteUrl "${PREFIX_CTEST_DROP_SITE}${endCDashUrl}")
  set(${PREFIX_CDASH_SITE_URL_OUT} "https://${cdashSiteUrl}" PARENT_SCOPE)
endfunction()


# @FUNCTION: tribits_get_cdash_build_url_from_parts()
#
# Create CDash index.php URL from the build parts.
#
# Usage::
#
#   tribits_get_cdash_build_url_from_parts(
#     INDEX_PHP_URL <indexPhpUrl>
#     PROJECT_NAME <projectName>
#     SITE_NAME <siteName>
#     BUILD_NAME <buildName>
#     BUILD_STAMP <buildStamp>
#     CDASH_BUILD_URL_OUT <cdashBuildUrlOut>
#     )
#
# Note that spaces are allowed ``<siteName>``, ``<buildName>`` or
# ``<buildStamp>`` and those will be handled correctly to produce a valid URL.
#
function(tribits_get_cdash_build_url_from_parts)
  # Get arguments
  cmake_parse_arguments(
    PREFIX #prefix
    "" #options
    "INDEX_PHP_URL;PROJECT_NAME;SITE_NAME;BUILD_NAME;BUILD_STAMP;CDASH_BUILD_URL_OUT" #one_value_keywords
    "" #multi_value_keywords
    ${ARGN}
    )
  # Do replacements for spaces and special chars in data
  tribits_replace_chars_for_url("${PREFIX_PROJECT_NAME}" project)
  tribits_replace_chars_for_url("${PREFIX_SITE_NAME}" site)
  tribits_replace_chars_for_url("${PREFIX_BUILD_NAME}" buildname)
  tribits_replace_chars_for_url("${PREFIX_BUILD_STAMP}" buildstamp)
  # Build the URL
  set(cdashIndexProj "${PREFIX_INDEX_PHP_URL}?project=${project}")
  set(filtersPreTxt "filtercount=3&showfilters=1&filtercombine=and")
  set(siteFlt "field1=site&compare1=61&value1=${site}")
  set(buildnameFlt "field2=buildname&compare2=61&value2=${buildname}")
  set(buildStampFlt "field3=buildstamp&compare3=61&value3=${buildstamp}")
  set(cdashBuildUrl
    "${cdashIndexProj}&${filtersPreTxt}&${siteFlt}&${buildnameFlt}&${buildStampFlt}")
  set(${PREFIX_CDASH_BUILD_URL_OUT} "${cdashBuildUrl}" PARENT_SCOPE)
endfunction()


# Replace chars in a regular string for usage in a URL with CDash
#
function(tribits_replace_chars_for_url  inputStr  outputStrForUrlOutVar)
  set(outputStrForUrl "${inputStr}")
  string(REPLACE " " "%20" outputStrForUrl "${outputStrForUrl}")
  string(REPLACE "+" "%2B" outputStrForUrl "${outputStrForUrl}")
  set(${outputStrForUrlOutVar} "${outputStrForUrl}" PARENT_SCOPE)
endfunction()

#  LocalWords:  GitHub tribits url buildname buildstamp
#  LocalWords:  tribits TRIBITS
#  LocalWords:  cmake CMake CMAKE
#  LocalWords:  ctest CTEST cdash CDash CDASH
#  LocalWords:  SUBSTRING
#  LocalWords:  endif  endfunction
