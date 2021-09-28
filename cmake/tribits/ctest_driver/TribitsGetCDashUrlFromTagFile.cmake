include(TribitsReadTagFile)


# @FUNCTION: tribits_get_build_url_and_write_to_file()
#
# Construct the build URL on CDash given the site name, buildname, and
# buildstamp (take from the TAG file) from inside of a running ctest -S
# program and optionally write it to a file as well.
#
# Usage::
#
#   tribits_get_build_url_and_write_to_file(
#     <cdashBuildUrlOut> [ <cdashBuldUrlFile> ] )
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
# If the file name argument ``<cdashBuldUrlFile>`` is non-empty, then that
# CDash URL will be written to the file as a single line.
#
function(tribits_get_build_url_and_write_to_file  cdashBuildUrlOut  cdashBuildUrlFile)
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
  set(cdashBuildUrl "https://${cdashBuildUrl}")
  if (cdashBuildUrlFile)
    file(WRITE "${cdashBuildUrlFile}" "${cdashBuildUrl}")
  endif()
  set(${cdashBuildUrlOut} "${cdashBuildUrl}" PARENT_SCOPE)
endfunction()


# @FUNCTION: tribits_print_build_url()
#
# Print the URL on CDash where build results can be found.
#
# Usage::
#
#   tribits_print_build_url( <msg> <cdashBuildUrl> )
#
function(tribits_print_build_url  msg  cdashBuildUrl)
  message("\n${msg}\n")
  message("    ${cdashBuildUrl}\n")
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
  string(FIND "${PREFIX_CTEST_DROP_LOCATION}" "?" endOfSubmitPhpIdx)
  string(SUBSTRING "${PREFIX_CTEST_DROP_LOCATION}" 0 ${endOfSubmitPhpIdx} submitPhpPart)
  string(REPLACE "submit.php" "index.php" indexPhpPart "${submitPhpPart}")
  set(indexPhpUrl "${PREFIX_CTEST_DROP_SITE}${indexPhpPart}")
  SET(${PREFIX_INDEX_PHP_URL_OUT} "${indexPhpUrl}" PARENT_SCOPE)
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
