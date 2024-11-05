
# Define the include dirs for the TESTONLY lib.
tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#
# Create a test-only b_mixed_lang lib which uses the optional test dependent
# MixedLang package
#

assert_defined(${PACKAGE_NAME}_ENABLE_MixedLang)
if (${PACKAGE_NAME}_ENABLE_MixedLang)

  include(${CMAKE_CURRENT_LIST_DIR}/ShowLibErrors.cmake)
  
  # Define the TESTONLY library
  tribits_add_library( b_mixed_lang
     TESTONLY
     HEADERS b_mixed_lang.hpp
     SOURCES b_mixed_lang.cpp
     DEPLIBS pws_b
     ${EXTRA_TAL_ARGS} # Used for testing and to demonstrate errors
     )

endif()

#
# Create a test-only lib that does not depend on any of the package's main
# libs.
#

tribits_add_library( b_test_utils
   TESTONLY
   HEADERS b_test_utils.hpp
   SOURCES b_test_utils.cpp
   )

