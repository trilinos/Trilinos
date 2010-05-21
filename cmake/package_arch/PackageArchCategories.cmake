INCLUDE(FindListElement)
INCLUDE(MessageWrapper)
INCLUDE(Join)


# Define the valid categories that will be recognized in the CATEGORIES keyword
SET(${PROJECT_NAME}_VALID_CATEGORIES  BASIC NIGHTLY PERFORMANCE)

# TODO: ABove, We may want only the final project to define these categories
# and not just be general categories for all projects based on ProjectArch.
# Given the logic below, only the categories BASIC and NIGHTLY are specially
# recognized.

# This is a string used in help and error messages
JOIN(${PROJECT_NAME}_VALID_CATEGORIES_STR ", "  FALSE ${${PROJECT_NAME}_VALID_CATEGORIES})

#
# Check for invalid categories called as:
#
#  PACKAGE_ARCH_GET_INVALID_CATEGORIES(INVLAID_CATEGORIES CATEGORIES_LIST)
#
#  The list of categories to check comes in through the ARGN list.
#
FUNCTION(PACKAGE_ARCH_GET_INVALID_CATEGORIES  INVALID_CATEGORIES_OUT)
  #MESSAGE("PACKAGE_ARCH_GET_INVALID_CATEGORIES: ${INVALID_CATEGORIES_OUT} ${ARGN}")
  SET(INVALID_CATEGORIES "")
  FOREACH(CATEGORY_IN ${ARGN})
    #PRINT_VAR(CATEGORY_IN)
    SET(FOUND_CATEGORY FALSE)
    FIND_LIST_ELEMENT(${PROJECT_NAME}_VALID_CATEGORIES ${CATEGORY_IN} FOUND_CATEGORY)
    IF (NOT FOUND_CATEGORY)
      #MESSAGE(STATUS "Not found in list of valid categories!")
      SET(INVALID_CATEGORIES ${INVALID_CATEGORIES} ${CATEGORY_IN})
    ENDIF()
    #PRINT_VAR(INVALID_CATEGORIES)
  ENDFOREACH()
  SET(${INVALID_CATEGORIES_OUT} ${INVALID_CATEGORIES} PARENT_SCOPE)
  #PRINT_VAR(${INVALID_CATEGORIES_OUT})
ENDFUNCTION()


#
# Assert there are no invalid categories called as:
#
#  PACKAGE_ARCH_ASSERT_VALID_CATEGORIES(CATEGORIES_LIST)
#
#  The list of categories to check comes in through the ARGN list.
#
FUNCTION(PACKAGE_ARCH_ASSERT_VALID_CATEGORIES)
  #MESSAGE("PACKAGE_ARCH_ASSERT_VALID_CATEGORIES: ${ARGN}")
  SET(INVALID_CATEGORIES "DUMMYCAT")
  PACKAGE_ARCH_GET_INVALID_CATEGORIES(INVALID_CATEGORIES ${ARGN})
  #PRINT_VAR(INVALID_CATEGORIES)
  IF (INVALID_CATEGORIES)
    MESSAGE_WRAPPER(SEND_ERROR "Error: The categories '${INVALID_CATEGORIES}' are not"
      " in the list of valid categories '${${PROJECT_NAME}_VALID_CATEGORIES_STR}'!")
  ENDIF()
ENDFUNCTION()
