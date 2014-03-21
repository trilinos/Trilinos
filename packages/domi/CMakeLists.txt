# -*- cmake -*-
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# Define the package
TRIBITS_PACKAGE(Domi DISABLE_STRONG_WARNINGS)

#
# Set up package-specific options
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_ABC
  HAVE_DOMI_ARRAY_BOUNDSCHECK
  "Enable Domi feature for array bounds checking."
  ${${PROJECT_NAME}_ENABLE_DEBUG}
  )

#
# Add the libraries, tests, and examples
ADD_SUBDIRECTORY(src)
TRIBITS_ADD_TEST_DIRECTORIES(test)

#
# Add the Parameter List documentation
TRIBITS_ADD_TEST_DIRECTORIES(doc/parameterList)

#
# Do standard postprocessing
TRIBITS_PACKAGE_POSTPROCESS()
