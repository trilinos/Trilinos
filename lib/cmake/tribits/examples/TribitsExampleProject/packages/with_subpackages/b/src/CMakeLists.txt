tribits_configure_file(${PACKAGE_NAME}_config.h)
tribits_include_directories(${CMAKE_CURRENT_BINARY_DIR})

include(${CMAKE_CURRENT_LIST_DIR}/ShowLibErrors.cmake)

tribits_include_directories(${CMAKE_CURRENT_BINARY_DIR})
tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})
tribits_add_library(pws_b
  SOURCES B.cpp
  HEADERS B.hpp ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  ADDED_LIB_TARGET_NAME_OUT pws_b_TARGET_NAME
  ${EXTRA_TAL_ARGS}
  )

tribits_verbose_print_var(pws_b_TARGET_NAME)

# Example of how to set properties on created target
set_target_properties( ${pws_b_TARGET_NAME} PROPERTIES LINKER_LANGUAGE CXX )
# NOTE: You don't have to set the LINKER_LANGUAGE
