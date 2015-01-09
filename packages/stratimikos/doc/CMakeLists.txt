# -*- cmake -*-

FILE( COPY
      ${CMAKE_CURRENT_SOURCE_DIR}/parameterList
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
    )

#
# Add the Parameter List documentation
TRIBITS_ADD_TEST_DIRECTORIES(parameterList)
