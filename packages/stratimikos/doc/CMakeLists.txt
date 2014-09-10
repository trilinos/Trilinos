
string(TOLOWER ${PACKAGE_NAME} LOWER_PACKAGE_NAME)

FILE( COPY 
      ${Trilinos_SOURCE_DIR}/packages/${LOWER_PACKAGE_NAME}/doc/parameterList
      DESTINATION ${Trilinos_BINARY_DIR}/packages/${PACKAGE_NAME}/doc)
#
# Add the Parameter List documentation
TRIBITS_ADD_TEST_DIRECTORIES(parameterList)
