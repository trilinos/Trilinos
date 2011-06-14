INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE( Trios ENABLE_SHADOWING_WARNINGS )

PACKAGE_ADD_DEBUG_OPTION()

PACKAGE_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

#
# B) Set up package-specific options
#

#
# C) Add the libraries, tests, and examples
#

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
