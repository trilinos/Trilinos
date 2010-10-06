INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)
INCLUDE(CombinedOption)

#
# A) Define the package
#

PACKAGE(MueLu)

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  MUELU_DEBUG
  "Build with debug code enabled."
  OFF )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(gallery)

#JJH this doesn't exist yet
#PACKAGE_ADD_TEST_DIRECTORIES(test)

#JJH this exists but isn't active yet
#PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
