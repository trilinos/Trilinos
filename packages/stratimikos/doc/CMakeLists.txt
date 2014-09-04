FILE( COPY ${Trilinos_SOURCE_DIR}/packages/stratimikos/doc/parameterList
      DESTINATION ${Trilinos_BINARY_DIR}/packages/stratimikos/doc)
#
# Add the Parameter List documentation
TRIBITS_ADD_TEST_DIRECTORIES(parameterList)
