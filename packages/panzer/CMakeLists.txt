INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)
INCLUDE(CombinedOption)

#
# A) Define the package
#

PACKAGE(Panzer)

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  PANZER_DEBUG
  "Build with debug code enabled."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  PANZER_STK_DEBUG
  "Build with debug code enabled."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Trios
  HAVE_IOSS
  "Build with IOSS code enabled."
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Trios
  ALBANY_IOSS
  "Build with Albany's IOSS code enabled."
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_EXPLICIT_TEMPLATE_INSTANTIATION
  "PANZER_EXPLICIT_TEMPLATE_INSTANTIATION"
  "Build with explicit template instantion enabled.  Defaults to ON."
  ON )

## Explicit instatiation support
#ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION
#  PANZER_HAVE_EXPLICIT_INSTANTIATION
#  "Build Panzer with explicit instantiation enabled."
#  OFF )

#ADD_OPTION(${PACKAGE_NAME}_EXPLICIT_INSTANTIATION_EVALUATION_TYPES
#  "Set the explicit instantiation evaluation types."
#  "panzer::Traits::Residual 2 panzer::Traits::Jacobian 1" )

#FILE(WRITE ${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE} ${DEPS_XML} )


#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(adapters)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
