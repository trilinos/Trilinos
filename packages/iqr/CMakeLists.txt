INCLUDE(TribitsPackageMacros)

#
# A) Forward declare the package so that certain options are also defined for
# subpackages
#

TRIBITS_PACKAGE_DECL(IQR) # ENABLE_SHADOWING_WARNINGS)

#
# B) Define the common options for iQR first so they can be used by
# subpackages as well.
#

TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_PACKAGE_DEF()

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)


TRIBITS_PACKAGE_POSTPROCESS()
