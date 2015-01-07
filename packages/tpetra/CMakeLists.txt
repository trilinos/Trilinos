############################################################
# Define the package
############################################################

# NOTE: If the TBB (Intel's Threading Building Blocks) TPL is enabled,
# it's probably not possible to build Tpetra with shadowing warnings.
#TRIBITS_PACKAGE( Tpetra ENABLE_SHADOWING_WARNINGS CLEANED )
TRIBITS_PACKAGE( Tpetra CLEANED )

TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

TRIBITS_SET_ST_FOR_DEV_MODE(ST_FOR_DEV_PT_FOR_RELEASE)

############################################################
# Process subpackages
############################################################

TRIBITS_PROCESS_SUBPACKAGES()

TRIBITS_PACKAGE_DEF()

#
# Identify what files to exclude when making a source tarball.
#
TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_EXCLUDE_FILES(
  core/doc/CodingGuidelines
  core/doc/TpetraDesign
  )

#
# Do standard package postprocessing
#
TRIBITS_PACKAGE_POSTPROCESS()
