INCLUDE(TribitsETISupport)

#
# Define the subpackage
#
TRIBITS_SUBPACKAGE(Kernels)

# Define ${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION,
# which defaults to ${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION.
TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

# Add this subpackage to the list of subpackages that do explicit
# template instantiation (ETI).
TRIBITS_ADD_ETI_SUPPORT()

# If building in debug mode, define the HAVE_TPETRAKERNELS_DEBUG macro.
TRIBITS_ADD_DEBUG_OPTION()

#
# Set up subpackage-specific configuration options
#

#
# "Optimization level" for TpetraKernels computational kernels.  The
# higher the level, the more code variants get generated, and thus the
# longer the compile times.  However, more code variants mean both
# better performance overall, and more uniform performance for corner
# cases.  Values of current interest (24 Apr 2014) are 0, 1, and 2.
#
TRIBITS_ADD_OPTION_AND_DEFINE( KokkosLinAlg_Opt_Level
  KOKKOSLINALG_OPT_LEVEL
  "Optimization level for TpetraKernels computational kernels: a nonnegative integer.  Higher levels result in better performance that is more uniform for corner cases, but increase build time and library size.  The default value is 1, which should give performance within ten percent of optimal on most platforms, for most problems."
  "1"
  )

# Enable experimental features of Tpetrakernels if set at configure
# time. Default is no.
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Experimental
  HAVE_TPETRAKERNELS_EXPERIMENTAL
  "Enable building and installation of experimental TpetraKernels features."
  NO )

# ==================================================================
# Explicit template instantiation (ETI) and test instantiation logic
# ==================================================================

# tpetra/CMakeLists.txt (the package's CMake logic) defines these
# variables.  Despite "ETI" in their names, they exist whether or not
# ETI is defined.  If ETI is defined, these variables govern the set
# of template parameter combinations over which Tpetra instantiates.
# Whether or not ETI is defined, the variables govern the set of
# template parameter combinations over which Tpetra runs tests.

ASSERT_DEFINED (Tpetra_ETI_SCALARS)
ASSERT_DEFINED (Tpetra_ETI_SCALARS_NO_ORDS)
ASSERT_DEFINED (Tpetra_ETI_LORDS)
ASSERT_DEFINED (Tpetra_ETI_DEVICES)

SET(${PACKAGE_NAME}_ETI_SCALARS "${Tpetra_ETI_SCALARS}")
# Exclude all ordinal types (GlobalOrdinal and int).
SET(${PACKAGE_NAME}_ETI_SCALARS_NO_ORDS "${Tpetra_ETI_SCALARS_NO_ORDS}")
SET(${PACKAGE_NAME}_ETI_LORDS "${Tpetra_ETI_LORDS}")
SET(${PACKAGE_NAME}_ETI_DEVICES "${Tpetra_ETI_DEVICES}")

# "Export" the names for use in the ETI system.
# If we don't do this, ETI won't see these variables.

GLOBAL_SET(${PACKAGE_NAME}_ETI_SCALARS ${${PACKAGE_NAME}_ETI_SCALARS})
GLOBAL_SET(${PACKAGE_NAME}_ETI_SCALARS_NO_ORDS ${${PACKAGE_NAME}_ETI_SCALARS_NO_ORDS})
GLOBAL_SET(${PACKAGE_NAME}_ETI_GORDS   ${${PACKAGE_NAME}_ETI_GORDS})
GLOBAL_SET(${PACKAGE_NAME}_ETI_LORDS   ${${PACKAGE_NAME}_ETI_LORDS})
GLOBAL_SET(${PACKAGE_NAME}_ETI_DEVICES ${${PACKAGE_NAME}_ETI_DEVICES})

# ==================================================================
# Process subdirectories
# ==================================================================

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(perf_test)
TRIBITS_ADD_TEST_DIRECTORIES(unit_test)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

TRIBITS_SUBPACKAGE_POSTPROCESS()

