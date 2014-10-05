#
# Create a test-only b_mixed_lang lib which uses the optional test dependent
# MixedLang package
#

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_MixedLanguage)
IF (${PACKAGE_NAME}_ENABLE_MixedLanguage)

  INCLUDE(${CMAKE_CURRENT_LIST_DIR}/ShowLibErrors.cmake)

  # Define the include dirs for the TESTONLY lib.
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  
  # Define the TESTONLY library
  TRIBITS_ADD_LIBRARY( b_mixed_lang
     TESTONLY
     HEADERS b_mixed_lang.hpp
     SOURCES b_mixed_lang.cpp
     DEPLIBS pws_b
     ${EXTRA_TAL_ARGS} # Used for testing and to demonstrate errors
     )

ENDIF()
