#
# A) Forward declare the package so that certain options are also defined for
# subpackages
#
TRIBITS_PACKAGE_DECL(WithSubpackages)

#
# B) Define the common options for the package first so they can be used by
# subpackages as well.
#
TRIBITS_ADD_DEBUG_OPTION()

#
# C) Process the subpackages
#
TRIBITS_PROCESS_SUBPACKAGES()

#
# D) Define the package now and perform standard post processing
#
TRIBITS_PACKAGE_DEF()
TRIBITS_PACKAGE_POSTPROCESS()
