
INCLUDE(ParseVariableArguments)
INCLUDE(PrintVar)


#
# Function that creates a set of tempalte code client headers that either
# includes or does not include the template instantiations depending on
# whether implicit or explicit instantiation is supported or not.
#
# Usage:
#
#    PACKAGE_CREATE_CLIENT_TEMPLATE_HEADERS(
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

FUNCTION(PACKAGE_CREATE_CLIENT_TEMPLATE_HEADERS BASE_DIR)

  #PRINT_VAR(BASE_DIR)
   
  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS(
    #prefix
    PARSE
    #lists
    ""
    #options
    "NOSIERRABJAM"
    ${ARGN}
    )

  #
  # B) Glob the names of all the X_decl.hpp files
  #

  FILE(GLOB DECL_HEADERS_LIST "${BASE_DIR}/*_decl.hpp")
  #PRINT_VAR(DECL_HEADERS_LIST)

  #
  # C) Write the client header files for each globed decl file
  #

  FOREACH(DECL_HEADER ${DECL_HEADERS_LIST})
 
    # Get the base file names (without _decl.hpp)
    STRING(REGEX REPLACE ".*/(.+)_decl.hpp" "\\1"  DECL_HEADER_BASE ${DECL_HEADER})
    #PRINT_VAR(DECL_HEADER_BASE)

    # Create the client header file
    SET(CLIENT_HEADER_STR "")
    APPEND_STRING_VAR(CLIENT_HEADER_STR
      "#include \"${DECL_HEADER_BASE}_decl.hpp\"\n"
       )
    IF (NOT HAVE_${PACKAGE_NAME_UC}_EXPLICIT_INSTANTIATION)
      APPEND_STRING_VAR(CLIENT_HEADER_STR
        "#include \"${DECL_HEADER_BASE}_def.hpp\"\n"
         )
    ENDIF()
    #PRINT_VAR(CLIENT_HEADER_STR)
    SET(BIN_HEADER_FILE "${CMAKE_CURRENT_BINARY_DIR}/${DECL_HEADER_BASE}.hpp")
    SET(WRITE_NEW_HEADER_FILE TRUE)
    IF (EXISTS "${BIN_HEADER_FILE}")
      # See if the file is the same and if it is, skip writing it again to avoid
      # unecessarily rebuilding object code.
      FILE(READ "${BIN_HEADER_FILE}" EXISTING_BIN_HEADER_STR)
      IF (CLIENT_HEADER_STR STREQUAL EXISTING_BIN_HEADER_STR)
        SET(WRITE_NEW_HEADER_FILE FALSE)
      ENDIF()
    ENDIF()
    IF (WRITE_NEW_HEADER_FILE)
      FILE(WRITE "${BIN_HEADER_FILE}" "${CLIENT_HEADER_STR}")
    ENDIF()

    # Create the SIERRA BJAM version of the header file
    IF (NOT PARSE_NOSIERRABJAM)
      SET(BJAM_CLIENT_HEADER_STR "")
      APPEND_STRING_VAR(BJAM_CLIENT_HEADER_STR
        "#include \"${DECL_HEADER_BASE}_decl.hpp\"\n"
        "#ifndef HAVE_${PACKAGE_NAME_UC}_EXPLICIT_INSTANTIATION\n"
        "#  include \"${DECL_HEADER_BASE}_def.hpp\"\n"
        "#endif\n"
         )
      SET(SIERRA_HEADER
        "${${PROJECT_NAME}_SOURCE_DIR}/SIERRA/bjam/config_headers/${DECL_HEADER_BASE}.hpp")
      IF (NOT EXISTS "${SIERRA_HEADER}")
        FILE(WRITE "${SIERRA_HEADER}" "${BJAM_CLIENT_HEADER_STR}")
      ENDIF()
    ENDIF()

  ENDFOREACH()

ENDFUNCTION()

