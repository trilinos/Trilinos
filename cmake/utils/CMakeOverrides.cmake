INCLUDE(ParseVariableArguments)

#This function is to override the standard behavior of include_directories.
#  We are overriding the default behavior for installation testing, this allows
#us to ensure that include directories will not be inadvertently added to the
#build lines for tests during installation testing. Normally we want the include
#directories to be handled as cmake usually does.However during installation
#testing we do not want most of the include directories to be used as the majority
#of the files should come from the installation we are building against. There is
#an exception to this and that is when there are test only headers that are needed.
#For that case we allow people to set "REQUIRED_DURING_INSTALLATION_TESTING" to
#tell us that this include directory does need to be set for instaltion testing.
FUNCTION(INCLUDE_DIRECTORIES)
  PARSE_ARGUMENTS(
    PARSE #prefix
    "" # Lists
    "REQUIRED_DURING_INSTALLATION_TESTING" #Options
    ${ARGN}
    )

  IF(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING OR PARSE_REQUIRED_DURING_INSTALLATION_TESTING)
    _INCLUDE_DIRECTORIES(${PARSE_DEFAULT_ARGS})
  ENDIF()
ENDFUNCTION()
