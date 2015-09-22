
INCLUDE(TribitsPackageMacros)

#
# A) Forward delcare the package so that certain options are also defined for
# subpackages
#

TRIBITS_PACKAGE_DECL(Pike ENABLE_SHADOWING_WARNINGS)

#
# B) Define the common options for Thyra first so they can be used by
# subpackages as well.
#

TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_ADD_ENABLE_TEUCHOS_TIME_MONITOR_OPTION()

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_EXPERIMENTAL
  HAVE_PIKE_EXPERIMENTAL
  "Enable experimental code not yet ready for release"
  OFF
  )

#
# C) Process the subpackages
#

TRIBITS_PROCESS_SUBPACKAGES()

#
# D) Process the package
#

TRIBITS_PACKAGE_DEF()

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_PACKAGE_POSTPROCESS()
