INCLUDE(Multiline_Set)
INCLUDE(Global_Set)
INCLUDE(Assert_Defined)
INCLUDE(Parse_Variable_Arguments)


MACRO(TPL_DECLARE_WITH_HEADER_LIBRARIES TPL_NAME)

  # Make sure the right name is used
  ASSERT_DEFINED(TPL_ENABLE_${TPL_NAME})

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "HEADER_NAMES;COMMON_LIB_NAMES"
     #options
     ""
     ${ARGN}
     )

  IF (TPL_${TPL_NAME}_LIBRARIES)
  
    # The ${TPL_NAME} library has already been found or the user has specified
    # these manually.  In this case, just make sure that everything has been
    # specified correctly
 
    ASSERT_DEFINED(TPL_${TPL_NAME}_INCLUDE_DIRS)

    IF (NOT DEFINED TPL_${TPL_NAME}_LIBRARY_DIRS)
      ADVANCED_SET(TPL_${TPL_NAME}_LIBRARY_DIRS)
    ENDIF()
 
    # Verify that indeed we have found the ${TPL_NAME} libraries!
  
    # ToDo: Implement!
  
  ELSE()
  
    # Otherwise, we need to look for ${TPL_NAME} library

    # Include directories

    IF (PARSE_HEADER_NAMES)
      SET(NAMES_ARG NAMES ${PARSE_HEADER_NAMES})
    ELSE()
      SET(NAMES_ARG)
    ENDIF()
  
    FIND_PATH(TPL_${TPL_NAME}_INCLUDE_DIRS NAMES ${NAMES_ARG}
      DOC "Full path to the ${TPL_NAME} headers")
    MARK_AS_ADVANCED(TPL_${TPL_NAME}_LIBRARIES)
  
    IF (NOT TPL_${TPL_NAME}_INCLUDE_DIRS)
      MULTILINE_SET(ERRMSG
        "Error, could not find the ${TPL_NAME} headers!  Please manually set"
        " TPL_${TPL_NAME}_INCLUDE_DIRS to point to the ${TPL_NAME} headers!")
      MESSAGE(FATAL_ERROR ${ERRMSG})
    ENDIF()

    # Libraries

    IF (PARSE_COMMON_LIB_NAMES)
      SET(NAMES_ARG NAMES ${PARSE_COMMON_LIB_NAMES})
    ELSE()
      SET(NAMES_ARG)
    ENDIF()
  
    FIND_LIBRARY(TPL_${TPL_NAME}_LIBRARIES ${NAMES_ARG}
      DOC "Full path to the ${TPL_NAME} library(s)")
    MARK_AS_ADVANCED(TPL_${TPL_NAME}_LIBRARIES)
  
    IF (NOT TPL_${TPL_NAME}_LIBRARIES)
      MULTILINE_SET(ERRMSG
        "Error, could not find the ${TPL_NAME} Library!  Please manually set"
        " TPL_${TPL_NAME}_LIBRARIES to point to the ${TPL_NAME} libraries!")
      MESSAGE(FATAL_ERROR ${ERRMSG})
    ENDIF()
  
  ENDIF()

  GLOBAL_NULL_SET(TPL_${TPL_NAME}_LIBRARY_DIRS)

ENDMACRO()
