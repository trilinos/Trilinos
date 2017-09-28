INCLUDE(TribitsCreateClientTemplateHeaders)
INCLUDE(CMakeParseArguments)

#
# Function that creates a set of template code client headers that either
# includes or does not include the template instantiations depending on
# whether implicit or explicit instantiation is supported or not.
#
# Usage:
#
#    TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(
#      BASE_DIR
#      [NOSIERRABJAM]
#      )
#
# The arguments are:
#
#  BASE_DIR
#
#    The base directory where files with the extension *_decl.hpp will be
#    globed for.
#
#  NOSIERRABJAM
#
#    If set, then the files will not be written in the SIERRA/bjam directory
#    for SIERRA to pick up.  This option would be used for template code that
#    is only used in tests and not in libraries.
#

FUNCTION(TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS  BASE_DIR)

  CMAKE_PARSE_ARGUMENTS(
    #prefix
    PARSE
    #options
    "NOSIERRABJAM"
    #one_value_keywords
    ""
    #multi_value_keywords
    ""
    ${ARGN}
    )

  IF (NOT PARSE_NOSIERRABJAM)
    # NOTE: The SIERRA integration system no longer uses this functionality to
    # install configured headers into a SIERRA/ directory.  Therefore, I am
    # commenting this out.
    #SET(ADDITIONAL_OUTPUT_DIRS_ARG
    #  ADDITIONAL_OUTPUT_DIRS
    #  "${${PROJECT_NAME}_SOURCE_DIR}/SIERRA/bjam/config_headers")
  ELSE()
    SET(ADDITIONAL_OUTPUT_DIRS_ARG)
  ENDIF()

  TRIBITS_CREATE_CLIENT_TEMPLATE_HEADERS(${BASE_DIR} ${ADDITIONAL_OUTPUT_DIRS_ARG})

ENDFUNCTION()

