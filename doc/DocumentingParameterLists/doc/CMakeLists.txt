FILE( COPY ${Trilinos_SOURCE_DIR}/packages/${PACKAGE_NAME}/doc/parameterList
      DESTINATION ${Trilinos_BINARY_DIR}/packages/${PACKAGE_NAME}/doc)
#
# Add the Parameter List documentation
TRIBITS_ADD_TEST_DIRECTORIES(parameterList)
