INCLUDE(TribitsPackageMacros)

#
# A) Forward delcare the package so that certain options are also defined for
# subpackages
#

TRIBITS_PACKAGE_DECL(Kokkos) # ENABLE_SHADOWING_WARNINGS)

#
# B) Define the common options for Kokkos first so they can be used by
# subpackages as well.
#

TRIBITS_ADD_DEBUG_OPTION()

#
# C) Process the subpackages for Kokkos
#

TRIBITS_PROCESS_SUBPACKAGES()

#
# D) If Kokkos itself is enabled, process the Kokkos package
#

TRIBITS_PACKAGE_DEF()

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_PACKAGE_POSTPROCESS()

