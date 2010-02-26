
INCLUDE(ParseVariableArguments)
INCLUDE(PrintVar)


#

# Function that creates a set of tempalte code client headers that either
# includes or does not include the template instantiations depending on
# whether implicit or explicit instantiation is supported or not.

FUNCTION(PACKAGE_CREATE_CLIENT_TEMPLATE_HEADERS  BASE_DIR)

  #PRINT_VAR(BASE_DIR)

  # Get the names of all the X_decl.hpp files
  FILE(GLOB DECL_HEADERS_LIST "${BASE_DIR}/*_decl.hpp")
  #PRINT_VAR(DECL_HEADERS_LIST)

  # Write the client header files
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
    FILE(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${DECL_HEADER_BASE}.hpp"
     ${CLIENT_HEADER_STR})

    # Create the SIERRA BJAM version of the header file
    SET(BJAM_CLIENT_HEADER_STR "")
    APPEND_STRING_VAR(BJAM_CLIENT_HEADER_STR
      "#include \"${DECL_HEADER_BASE}_decl.hpp\"\n"
      "#ifndef HAVE_${PACKAGE_NAME_UC}_EXPLICIT_INSTANTIATION\n"
      "#  include \"${DECL_HEADER_BASE}_def.hpp\"\n"
      "#endif\n"
       )
    FILE(WRITE
     "${${PROJECT_NAME}_SOURCE_DIR}/SIERRA/bjam/config_headers/${DECL_HEADER_BASE}.hpp"
     ${BJAM_CLIENT_HEADER_STR})

  ENDFOREACH()

ENDFUNCTION()

